Yoshida_2018_coproporphyrin_I_rifampin <- function() {
  description <- "One-compartment endogenous turnover model for the OATP1B-substrate biomarker coproporphyrin I (CPI) in healthy adults (Yoshida 2018, rifampin-CPI calibration). CPI is produced at a zero-order synthesis rate Ksyn = kdeg * Baseline and eliminated as a single first-order pool whose overall rate constant kdeg is decomposed into a non-hepatic fraction fNH (unaffected by inhibitor) and a hepatic fraction 1 - fNH (competitively inhibited by the OATP1B perpetrator via Ki,u). The perpetrator portal-vein unbound concentration enters as a time-varying covariate CP_RIF_UM (umol/L); setting CP_RIF_UM = 0 collapses the model to the inhibitor-free steady state Baseline. This file encodes the rifampin-CPI calibration (Table 2 left column; no IIV reported); a sibling file Yoshida_2018_coproporphyrin_I_GDC0810 encodes the GDC-0810 calibration with its own Ki,u, kdeg, and IIV structure. The original fit used a Simcyp v16r1 default single-dose rifampin model for the portal-vein concentration profile; that PBPK output is not reproducible from on-disk sources and the paper itself documents an approximately 5-fold sensitivity of the estimated Ki,u to the choice of perpetrator-PK model, so downstream users must supply CP_RIF_UM externally and treat the calibrated Ki,u as conditional on that choice."
  reference <- paste(
    "Yoshida K, Guo C, Sane R.",
    "Quantitative Prediction of OATP-Mediated Drug-Drug Interactions",
    "With Model-Based Analysis of Endogenous Biomarker Kinetics.",
    "CPT Pharmacometrics Syst. Pharmacol. 2018;7(8):517-524.",
    "doi:10.1002/psp4.12315.",
    "The rifampin-CPI calibration (Table 2 left column) was fit",
    "to the Lai et al. 2016 plasma CPI profile cohort using",
    "portal-vein rifampin concentrations from the Simcyp v16r1",
    "default single-dose rifampin model file as the forcing",
    "function; that Simcyp output is not reproducible from on-disk",
    "sources, so users must supply CP_RIF_UM externally.",
    "The companion GDC-0810 calibration is parameterised in",
    "modellib('Yoshida_2018_coproporphyrin_I_GDC0810').",
    sep = " "
  )
  vignette <- "Yoshida_2018_coproporphyrin_I_rifampin"
  units <- list(time = "hour", dosing = "none", concentration = "nmol/L")

  covariateData <- list(
    CP_RIF_UM = list(
      description        = "Instantaneous rifampin portal-vein unbound concentration as a time-varying perpetrator covariate driving competitive OATP1B inhibition of the hepatic component of CPI clearance (Yoshida 2018 Methods, Model-based analysis with inhibitor kinetics).",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Set to 0 outside the rifampin dosing window so the hepatic-inhibition term collapses to the no-inhibition form and the model returns to the steady-state Baseline. The original Yoshida 2018 fit used the Simcyp v16r1 default single-dose rifampin model file output for portal-vein unbound concentration; that PBPK profile is not reproducible from on-disk sources. The paper itself reports (Discussion) that a preliminary refit with a different rifampin PK model that matched an alternative observed time-course produced an estimated Ki,u about 5-fold higher, so the calibrated Ki,u below is conditional on the choice of perpetrator-PK model. Users must supply CP_RIF_UM externally.",
      source_name        = "CRIF"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 12L,
    n_studies       = 1L,
    age_range       = "Healthy adult subjects; per-subject demographics not tabulated by Yoshida 2018. The underlying clinical dataset is Lai et al. 2016 (J Pharmacol Exp Ther 358:397-404), 12 healthy male SLCO1B1 c.521 T>C wildtype subjects.",
    weight_range    = "(not extracted; Yoshida 2018 does not tabulate per-subject weights for the rifampin-CPI cohort.)",
    sex_female_pct  = 0,
    disease_state   = "Healthy male adult volunteers in a single-dose rifampin drug-drug-interaction study; CPI monitored as a candidate endogenous biomarker of OATP1B-mediated DDIs.",
    dose_range      = "Endogenous biomarker (no exogenous CPI dose); rifampin co-administered as a 600 mg oral dose (per Lai et al. 2016 study design).",
    regions         = "(not extracted; Yoshida 2018 does not state the region for the Lai et al. 2016 cohort.)",
    notes           = "Demographics inferred from Yoshida 2018 Methods and the cited source clinical study Lai et al. 2016. Yoshida 2018 reports kinetic profiles of CPI 'in the presence of rifampin' as the input data for the model-based analysis; the rifampin-CPI fit (Table 2 left column) does not report IIV, consistent with a typical-value fit of the n=12 cohort."
  )

  ini({
    # Structural parameters -- Table 2, RIF-CPI interaction column.
    # The one-compartment CPI model has only four structural parameters
    # (Baseline, kdeg, fNH, Ki,u) plus a proportional residual error.
    lrbase   <- log(0.863)
    label("CPI baseline plasma concentration (nmol/L)")
    # Table 2, RIF-CPI: Baseline = 0.863 nM (RSE 4.61%); steady-state
    # plasma CPI with no inhibitor present.

    lkdeg    <- log(2.55)
    label("CPI overall first-order degradation rate constant (1/h)")
    # Table 2, RIF-CPI: kdeg = 2.55 hr^-1 (RSE 8.88%); the total
    # apparent elimination rate constant decomposed into hepatic and
    # non-hepatic components via fNH.

    logitfnh <- qlogis(0.129)
    label("Logit of the non-hepatic fraction fNH of overall CPI clearance (unitless; bounded to [0,1])")
    # Table 2, RIF-CPI: fNH = 12.9 % (RSE 6.66%), estimated; the
    # complementary 1 - fNH = 87.1 % is the OATP1B-mediated hepatic
    # fraction subject to perpetrator inhibition. Logit-transformed
    # in ini() to enforce 0 <= fNH <= 1 in any downstream re-fitting.

    lkiu     <- log(0.0203)
    label("Rifampin unbound OATP1B inhibition constant Ki,u (umol/L)")
    # Table 2, RIF-CPI: Ki,u = 0.0203 uM (RSE 17.0%); conditional on
    # the Simcyp v16r1 default rifampin portal-vein concentration
    # profile used as the forcing function. Yoshida 2018 Discussion
    # documents an approximately 5x sensitivity of this estimate to
    # the choice of rifampin PBPK model.

    # Residual error -- Table 2 reports a proportional residual error.
    propSd   <- 0.0513
    label("Proportional residual error (fraction)")
    # Table 2, RIF-CPI: proportional residual error = 5.13 %CV
    # (RSE 20.4 %) -> fraction 0.0513.
  })

  model({
    # Individual parameters (no IIV reported in Table 2 for the
    # RIF-CPI fit; all etas omitted intentionally).
    rbase <- exp(lrbase)
    kdeg  <- exp(lkdeg)
    fnh   <- expit(logitfnh)
    kiu   <- exp(lkiu)

    # Synthesis rate is anchored to the steady-state identity at no
    # inhibitor: Ksyn = kdeg * Baseline.
    ksyn <- kdeg * rbase

    # Competitive OATP1B inhibition reduces the hepatic component of
    # kdeg only; the non-hepatic fraction fNH is unaffected by the
    # perpetrator. CP_RIF_UM = 0 collapses kdeg_eff to kdeg.
    kdeg_eff <- kdeg * (fnh + (1 - fnh) / (1 + CP_RIF_UM / kiu))

    # ODE: one-compartment turnover. central is the plasma CPI
    # concentration in nmol/L; no volume of distribution appears
    # explicitly because the Yoshida 2018 parameterisation operates
    # directly on the concentration (Ksyn carries units of nmol/L/h
    # so dimensional consistency is satisfied).
    d/dt(central) <- ksyn - kdeg_eff * central
    central(0)    <- rbase

    Cc <- central
    Cc ~ prop(propSd)
  })
}
