Zhang_2025_remdesivir_esrd <- function() {
  description <- paste(
    "Six-state mechanism-based cascade model for intravenous remdesivir and",
    "its two circulating metabolites GS-704277 (the paper's 'IM') and",
    "GS-441524 (the paper's 'NUC'), refitted to a single anuric patient with",
    "end-stage renal disease (eGFR = 0 mL/min) on intermittent haemodialysis",
    "(Zhang 2025, 'Clinical study based simulation 2'). The structure is",
    "identical to the healthy-and-renally-impaired companion model",
    "Zhang_2025_remdesivir, but every parameter was re-estimated against a",
    "different data source -- the Soergel 2021 case report -- so the two are",
    "separate fits rather than two arms of one covariate model. Two",
    "differences are structural rather than numerical: central conversion of",
    "remdesivir to GS-704277 is estimated at exactly zero, so in this",
    "patient the whole GS-704277 formation flux is peripheral, and the",
    "GS-441524-to-triphosphate rate constant is an ordinary small number",
    "(0.037 /h) rather than the effectively infinite sink the cohort model",
    "carries. Written entirely in CONCENTRATION space: every state is a",
    "concentration in ng/mL and every parameter a first-order rate constant",
    "in 1/h, so no volume of distribution appears and a dose must be",
    "supplied as the initial central remdesivir concentration.")
  reference <- "Zhang S, Jeong S, Jiang B, Ho H. Pharmacokinetic simulations for remdesivir and its metabolites in healthy subjects and patients with renal impairment. Front Pharmacol. 2025;16:1488961. doi:10.3389/fphar.2025.1488961"
  vignette <- "Zhang_2025_remdesivir"

  # `dosing` is deliberately ng/mL, not mg -- see the companion model
  # Zhang_2025_remdesivir and the `compartmentData` note below.
  units <- list(time = "h", dosing = "ng/mL", concentration = "ng/mL")

  # Every state holds a CONCENTRATION (ng/mL), not an amount, exactly as in
  # Zhang 2025 Equations 1-6. No volume of distribution is reported
  # anywhere in the paper, so a dose enters as an initial central
  # remdesivir concentration.
  # Peripheral compartments recorded as `blood cell`: Zhang 2025 identifies
  # them with peripheral blood mononuclear cells (Introduction and Discussion).
  # See the companion model Zhang_2025_remdesivir for the full quotation.
  compartmentData <- list(
    central                = list(analyte = "remdesivir (GS-5734)", units = "ng/mL", specimen = "plasma", verified = TRUE),
    peripheral1            = list(analyte = "remdesivir (GS-5734)", units = "ng/mL", specimen = "blood cell", verified = TRUE),
    central_gs704277       = list(analyte = "GS-704277", units = "ng/mL", specimen = "plasma", verified = TRUE),
    peripheral1_gs704277   = list(analyte = "GS-704277", units = "ng/mL", specimen = "blood cell", verified = TRUE),
    central_gs441524       = list(analyte = "GS-441524", units = "ng/mL", specimen = "plasma", verified = TRUE),
    peripheral1_gs441524   = list(analyte = "GS-441524", units = "ng/mL", specimen = "blood cell", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 1L,
    n_studies      = 1L,
    age_range      = "A male patient in his mid-seventies (Zhang 2025 Methods 2.1; the source case report withholds the exact age).",
    sex_female_pct = 0,
    disease_state  = "A single kidney-transplant recipient in his mid-seventies with COVID-19, receiving renal replacement therapy and with an eGFR of 0 mL/min, i.e. no residual renal function. Described in the case report of Soergel F, Malin JJ, Hagmann H, et al., J Antimicrob Chemother. 2021;76:825-827, doi:10.1093/jac/dkaa500, and used by Zhang 2025 as an extreme-renal-impairment test of the model.",
    renal_function = "eGFR 0 mL/min (anuric), on intermittent haemodialysis. This is well below the eGFR 30 mL/min/1.73 m^2 threshold at which remdesivir is contraindicated in labelling.",
    dose_range     = "Standard 5-day regimen: a 200 mg intravenous infusion on day 1 followed by 100 mg daily on each of the next 4 days. Zhang 2025 fitted and validated against DAY 1 ONLY ('For model validation, we specifically analyzed data from the first day of administration'), so the parameters describe a single 200 mg dose and the model is uninformed about accumulation over the remaining four days or about the effect of the intervening dialysis sessions.",
    notes          = "IMPORTANT PROVENANCE LIMITATION. As for the companion cohort model, Zhang 2025 had no individual-level data: the concentration-time profiles were digitised from the published case-report figures with Engauge Digitizer 12.1 (Methods 2.1). Here that means the entire dataset is one patient's digitised curves for three analytes, so the estimates are a deterministic curve fit rather than a population analysis, and no variability of any kind is estimable. Zhang 2025 Discussion notes that no GS-443902 (triphosphate) data were reported for this patient either."
  )

  ini({
    # ------------------------------------------------------------------
    # STRUCTURAL PARAMETERS -- eGFR = 0 column of Zhang 2025 Table 1
    # ("Parameters for a renal-impaired patient with eGFR = 0
    # (Soergel et al., 2021)").
    #
    # As in the companion model, every Table 1 parameter carries units of
    # /h and is named here for the micro-constant it is rather than for the
    # clearance-flavoured symbol Table 1 uses. Point estimates only; no
    # standard errors or confidence intervals are reported.
    # ------------------------------------------------------------------

    lk12 <- log(0.31)
    label("Remdesivir central-peripheral exchange rate constant (1/h)")                      # Zhang 2025 Table 1 row Q_RDV, eGFR = 0 column: 0.31 /h

    lkel <- log(2.3)
    label("Remdesivir elimination rate constant from the central compartment (1/h)")         # Zhang 2025 Table 1 row CL_C,RDV, eGFR = 0 column: 2.3 /h. Table 1 calls this "Clearance of RDV" but its units are /h, so in this concentration-space formulation it is an elimination rate constant.

    lk12_gs704277 <- log(83.8)
    label("GS-704277 central-peripheral exchange rate constant (1/h)")                       # Zhang 2025 Table 1 row Q_IM, eGFR = 0 column: 83.8 /h

    lk12_gs441524 <- log(50.02)
    label("GS-441524 central-peripheral exchange rate constant (1/h)")                       # Zhang 2025 Table 1 row Q_NUC, eGFR = 0 column: 50.02 /h

    # Central remdesivir-to-GS-704277 conversion is estimated at EXACTLY
    # ZERO in this patient, so in the eGFR = 0 fit all GS-704277 is formed
    # in the peripheral compartment. This is the one parameter in the file
    # that is NOT log-transformed: log(0) is -Inf and cannot be carried as
    # a population parameter, so the value is held on the linear scale and
    # used directly in the ODE below. It is wrapped in `fixed()` because a
    # boundary estimate of exactly zero is a structural statement about the
    # model rather than an interior point estimate.
    kmet_gs704277_central <- fixed(0)
    label("Remdesivir-to-GS-704277 conversion rate constant, central compartment (1/h)")     # Zhang 2025 Table 1 row K_C,IM, eGFR = 0 column: 0 /h

    lkmet_gs704277_peripheral1 <- log(0.22)
    label("Remdesivir-to-GS-704277 conversion rate constant, peripheral compartment (1/h)")  # Zhang 2025 Table 1 row K_P,IM, eGFR = 0 column: 0.22 /h

    lkmet_gs441524_central <- log(6.27)
    label("GS-704277-to-GS-441524 conversion rate constant, central compartment (1/h)")      # Zhang 2025 Table 1 row K_C,NUC, eGFR = 0 column: 6.27 /h

    lkmet_gs441524_peripheral1 <- log(0.15)
    label("GS-704277-to-GS-441524 conversion rate constant, peripheral compartment (1/h)")   # Zhang 2025 Table 1 row K_P,NUC, eGFR = 0 column: 0.15 /h

    # Unlike the cohort model, where this rate constant is ~3e10 /h and
    # acts as an instantaneous sink, here it is an ordinary small number,
    # so the peripheral GS-441524 pool is genuinely retained and the ODE
    # system is not stiff.
    lkmet_gs443902_peripheral1 <- log(0.037)
    label("GS-441524-to-GS-443902 conversion rate constant, peripheral compartment (1/h)")   # Zhang 2025 Table 1 row K_P,NTP, eGFR = 0 column: 0.037 /h

    # ------------------------------------------------------------------
    # NO INTER-INDIVIDUAL VARIABILITY IS CARRIED, and none is estimable:
    # this fit has a single subject. Zhang 2025 Methods 2.3 declares
    # lognormal random effects on all parameters for the mixed-effects
    # analysis, but reports no variance, SD or CV% for any parameter in
    # either fit, and there is no supplement (the EuropePMC
    # supplementaryFiles bundle for PMC11982744 contains only the six
    # publisher figure files). Etas are OMITTED rather than written as
    # `~ fixed(0)` because a zero-variance diagonal makes OMEGA singular
    # and breaks the Cholesky sampler used by rxSolve. Recorded in the
    # vignette Errata.
    #
    # RESIDUAL ERROR: Zhang 2025 Equation 7 declares a combined
    # additive-plus-proportional model but reports neither a nor b, so
    # both terms are fixed at zero for all three analytes and simulations
    # from this file are noise-free.
    # ------------------------------------------------------------------

    addSd <- fixed(0)
    label("Additive residual SD on remdesivir (ng/mL); magnitude not reported")              # Zhang 2025 Equation 7 declares the combined error model but reports no value for a
    propSd <- fixed(0)
    label("Proportional residual SD on remdesivir (fraction); magnitude not reported")       # Zhang 2025 Equation 7 declares the combined error model but reports no value for b

    addSd_gs704277 <- fixed(0)
    label("Additive residual SD on GS-704277 (ng/mL); magnitude not reported")               # Zhang 2025 Equation 7 declares the combined error model but reports no value for a
    propSd_gs704277 <- fixed(0)
    label("Proportional residual SD on GS-704277 (fraction); magnitude not reported")        # Zhang 2025 Equation 7 declares the combined error model but reports no value for b

    addSd_gs441524 <- fixed(0)
    label("Additive residual SD on GS-441524 (ng/mL); magnitude not reported")               # Zhang 2025 Equation 7 declares the combined error model but reports no value for a
    propSd_gs441524 <- fixed(0)
    label("Proportional residual SD on GS-441524 (fraction); magnitude not reported")        # Zhang 2025 Equation 7 declares the combined error model but reports no value for b
  })

  model({
    # ------------------------------------------------------------------
    # Individual parameters. No covariates: this is a single-patient fit,
    # so every parameter is used at its typical value.
    # `kmet_gs704277_central` is already on the linear scale (fixed at 0)
    # and is used directly in equations 1 and 3 without back-transforming.
    # ------------------------------------------------------------------
    k12 <- exp(lk12)
    kel <- exp(lkel)
    k12_gs704277 <- exp(lk12_gs704277)
    k12_gs441524 <- exp(lk12_gs441524)
    kmet_gs704277_peripheral1 <- exp(lkmet_gs704277_peripheral1)
    kmet_gs441524_central <- exp(lkmet_gs441524_central)
    kmet_gs441524_peripheral1 <- exp(lkmet_gs441524_peripheral1)
    kmet_gs443902_peripheral1 <- exp(lkmet_gs443902_peripheral1)

    # ------------------------------------------------------------------
    # ODE system -- Zhang 2025 Equations 1-6, transcribed one for one.
    # Identical structure to the companion cohort model; only the
    # parameter values differ. See Zhang_2025_remdesivir for the full
    # equation listing and the discussion of the model's asymmetries
    # (remdesivir eliminated only centrally, GS-441524 removed only
    # peripherally).
    # ------------------------------------------------------------------
    d/dt(central) <-
      k12 * (peripheral1 - central) -
      kmet_gs704277_central * central -
      kel * central                                                          # Equation 1
    d/dt(peripheral1) <-
      k12 * (central - peripheral1) -
      kmet_gs704277_peripheral1 * peripheral1                                # Equation 2

    d/dt(central_gs704277) <-
      k12_gs704277 * (peripheral1_gs704277 - central_gs704277) +
      kmet_gs704277_central * central -
      kmet_gs441524_central * central_gs704277                               # Equation 3
    d/dt(peripheral1_gs704277) <-
      k12_gs704277 * (central_gs704277 - peripheral1_gs704277) +
      kmet_gs704277_peripheral1 * peripheral1 -
      kmet_gs441524_peripheral1 * peripheral1_gs704277                       # Equation 4

    d/dt(central_gs441524) <-
      k12_gs441524 * (peripheral1_gs441524 - central_gs441524) +
      kmet_gs441524_central * central_gs704277                               # Equation 5
    d/dt(peripheral1_gs441524) <-
      k12_gs441524 * (central_gs441524 - peripheral1_gs441524) +
      kmet_gs441524_peripheral1 * peripheral1_gs704277 -
      kmet_gs443902_peripheral1 * peripheral1_gs441524                       # Equation 6

    # ------------------------------------------------------------------
    # Observations. The states already hold plasma concentrations in
    # ng/mL, so no division by a volume is required. Zhang 2025 Figure 5
    # plots these three central-compartment quantities against the
    # digitised Soergel 2021 observations, with ng/mL on every y axis.
    # ------------------------------------------------------------------
    Cc <- central
    Cc_gs704277 <- central_gs704277
    Cc_gs441524 <- central_gs441524

    Cc ~ add(addSd) + prop(propSd)
    Cc_gs704277 ~ add(addSd_gs704277) + prop(propSd_gs704277)
    Cc_gs441524 ~ add(addSd_gs441524) + prop(propSd_gs441524)
  })
}
