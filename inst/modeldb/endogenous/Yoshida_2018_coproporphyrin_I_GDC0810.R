Yoshida_2018_coproporphyrin_I_GDC0810 <- function() {
  description <- "One-compartment endogenous turnover model for the OATP1B-substrate biomarker coproporphyrin I (CPI) in healthy adults (Yoshida 2018, GDC-0810-CPI calibration). CPI is produced at a zero-order synthesis rate Ksyn = kdeg * Baseline and eliminated as a single first-order pool whose overall rate constant kdeg is decomposed into a non-hepatic fraction fNH (held fixed at 12.9 %, unaffected by inhibitor) and a hepatic fraction 1 - fNH (competitively inhibited by the OATP1B perpetrator via Ki,u). The perpetrator portal-vein unbound concentration enters as a time-varying covariate CP_GDC_UM (umol/L); setting CP_GDC_UM = 0 collapses the model to the inhibitor-free steady state Baseline. This file encodes the GDC-0810-CPI calibration (Table 2 right column) with IIV on Baseline (18.2 %CV) and Ki,u (30.1 %CV); a sibling file Yoshida_2018_coproporphyrin_I_rifampin encodes the rifampin calibration with its own Ki,u, kdeg, and no IIV. The original fit used a Y. Chen et al. in-house PBPK model for GDC-0810 portal-vein concentrations (personal communication, not on disk and not in the nlmixr2lib registry), so downstream users must supply CP_GDC_UM externally."
  reference <- paste(
    "Yoshida K, Guo C, Sane R.",
    "Quantitative Prediction of OATP-Mediated Drug-Drug Interactions",
    "With Model-Based Analysis of Endogenous Biomarker Kinetics.",
    "CPT Pharmacometrics Syst. Pharmacol. 2018;7(8):517-524.",
    "doi:10.1002/psp4.12315.",
    "The GDC-0810-CPI calibration (Table 2 right column) was fit",
    "to the Liu et al. 2018 plasma CPI profile cohort using",
    "portal-vein GDC-0810 concentrations from an in-house Y. Chen",
    "et al. PBPK model (cited as personal communication; not on",
    "disk), so users must supply CP_GDC_UM externally. The",
    "companion rifampin calibration is parameterised in",
    "modellib('Yoshida_2018_coproporphyrin_I_rifampin').",
    sep = " "
  )
  vignette <- "Yoshida_2018_coproporphyrin_I_GDC0810"
  units <- list(time = "hour", dosing = "none", concentration = "nmol/L")

  covariateData <- list(
    CP_GDC_UM = list(
      description        = "Instantaneous GDC-0810 portal-vein unbound concentration as a time-varying perpetrator covariate driving competitive OATP1B inhibition of the hepatic component of CPI clearance (Yoshida 2018 Methods, Model-based analysis with inhibitor kinetics).",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Set to 0 outside the GDC-0810 dosing window so the hepatic-inhibition term collapses to the no-inhibition form and the model returns to the steady-state Baseline. The original Yoshida 2018 fit used an in-house Y. Chen et al. PBPK model output for GDC-0810 portal-vein unbound concentration (cited as personal communication in the paper); that profile is not reproducible from on-disk sources and no GDC-0810 PK model is currently registered in nlmixr2lib. Users must supply CP_GDC_UM externally. The paper notes (Discussion) that observed GDC-0810 plasma AUC IIV was about 20 %, so the IIV reported on Ki,u below (30.1 %CV) partially includes per-subject variability in portal-vein exposure rather than purely intrinsic Ki,u variability.",
      source_name        = "CGDC"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = NA_integer_,
    n_studies       = 1L,
    age_range       = "Healthy adult female subjects; per-subject demographics not tabulated by Yoshida 2018. The underlying clinical dataset is Liu et al. 2018 (J Clin Pharmacol; effect of OATP1B1/1B3 inhibitor GDC-0810 on the pharmacokinetics of pravastatin and coproporphyrin I/III in healthy female subjects).",
    weight_range    = "(not extracted; Yoshida 2018 does not tabulate per-subject weights for the GDC-0810-CPI cohort.)",
    sex_female_pct  = 100,
    disease_state   = "Healthy female adult volunteers in a GDC-0810 drug-drug-interaction study; CPI monitored as an endogenous biomarker of OATP1B-mediated DDIs.",
    dose_range      = "Endogenous biomarker (no exogenous CPI dose); GDC-0810 was administered orally per the Liu et al. 2018 study design.",
    regions         = "(not extracted; Yoshida 2018 does not state the region for the Liu et al. 2018 cohort.)",
    notes           = "Demographics inferred from Yoshida 2018 Methods and the cited source clinical study Liu et al. 2018. Yoshida 2018 reports that 'both the central trend and variability were described with the estimated parameters including interindividual variability (IIV)' for the GDC-0810-CPI fit, consistent with the IIV estimates on Baseline and Ki,u in Table 2 right column."
  )

  ini({
    # Structural parameters -- Table 2, GDC-0810-CPI interaction column.
    # Same one-compartment turnover structure as the rifampin sibling
    # (modellib('Yoshida_2018_coproporphyrin_I_rifampin')); parameter
    # values differ per calibration.
    lrbase   <- log(0.873)
    label("CPI baseline plasma concentration (nmol/L)")
    # Table 2, GDC-0810-CPI: Baseline = 0.873 nM (RSE 7.47%); steady-
    # state plasma CPI with no inhibitor present.

    lkdeg    <- log(1.25)
    label("CPI overall first-order degradation rate constant (1/h)")
    # Table 2, GDC-0810-CPI: kdeg = 1.25 hr^-1 (RSE 5.88%). Yoshida
    # 2018 Results notes this estimate is "comparable" to the
    # rifampin-CPI estimate (2.55 hr^-1) although the two values
    # differ approximately 2-fold; the authors interpret the kdeg
    # estimates as plausibly transferable across DDI scenarios.

    logitfnh <- fixed(qlogis(0.129))
    label("Logit of the non-hepatic fraction fNH of overall CPI clearance (unitless; fixed)")
    # Table 2, GDC-0810-CPI: fNH = 12.9 % FIXED (not estimated). The
    # authors fix fNH to the rifampin-CPI estimated value because
    # their sensitivity analysis showed fNH has small influence on
    # the other parameters in the GDC-0810-CPI fit. Wrapped in
    # fixed() to mark the structural-assumption provenance per
    # SKILL.md.

    lkiu     <- log(0.00174)
    label("GDC-0810 unbound OATP1B inhibition constant Ki,u (umol/L)")
    # Table 2, GDC-0810-CPI: Ki,u = 0.00174 uM (RSE 27.3%); about
    # 12x lower than the rifampin-CPI Ki,u, consistent with GDC-0810
    # being a more potent OATP1B inhibitor on an unbound-Ki basis.
    # Conditional on the Y. Chen et al. PBPK model for GDC-0810
    # portal-vein concentration (personal communication; not on disk).

    # IIV -- Table 2 reports IIV on Baseline and on Ki,u for the
    # GDC-0810-CPI fit; the rifampin-CPI fit has no IIV. Variances
    # computed via omega^2 = log(1 + CV^2).
    etalrbase ~ 0.03244
    # Table 2, GDC-0810-CPI: IIV on Baseline = 18.2 %CV (RSE 23.2%);
    # log(1 + 0.182^2) = 0.03244.

    etalkiu   ~ 0.08660
    # Table 2, GDC-0810-CPI: IIV on Ki,u = 30.1 %CV (RSE 54.3%);
    # log(1 + 0.301^2) = 0.08660. Per Yoshida 2018 Discussion this
    # likely includes per-subject variability in GDC-0810 portal-vein
    # exposure (about 20 %CV) in addition to true Ki,u variability.

    # Residual error -- Table 2 reports a proportional residual error.
    propSd   <- 0.119
    label("Proportional residual error (fraction)")
    # Table 2, GDC-0810-CPI: proportional residual error = 11.9 %CV
    # (RSE 14.9 %) -> fraction 0.119.
  })

  model({
    # Individual parameters with log-normal IIV on Baseline and Ki,u;
    # logit-scale fNH has no IIV (fixed structurally) and kdeg has no
    # IIV (no IIV reported for kdeg in Table 2).
    rbase <- exp(lrbase + etalrbase)
    kdeg  <- exp(lkdeg)
    fnh   <- expit(logitfnh)
    kiu   <- exp(lkiu + etalkiu)

    # Synthesis rate is anchored to the per-subject steady-state
    # identity at no inhibitor: Ksyn_i = kdeg * Baseline_i.
    ksyn <- kdeg * rbase

    # Competitive OATP1B inhibition reduces the hepatic component of
    # kdeg only; the non-hepatic fraction fNH is unaffected by the
    # perpetrator. CP_GDC_UM = 0 collapses kdeg_eff to kdeg.
    kdeg_eff <- kdeg * (fnh + (1 - fnh) / (1 + CP_GDC_UM / kiu))

    # ODE: one-compartment turnover. central is the plasma CPI
    # concentration in nmol/L; no volume of distribution appears
    # explicitly because the Yoshida 2018 parameterisation operates
    # directly on the concentration.
    d/dt(central) <- ksyn - kdeg_eff * central
    central(0)    <- rbase

    Cc <- central
    Cc ~ prop(propSd)
  })
}
