Ouerdani_2015_pazopanib <- function() {
  description <- "Semi-mechanistic tumour growth and angiogenesis-inhibition (TGI) model for pazopanib in renal-cell carcinoma patients (Ouerdani 2015 clinical fit): logistic tumour growth (state tumor_size) limited by a separately tracked vasculature-determined carrying capacity (state carrying_capacity), with antiangiogenic and cytotoxic drug effects parameterised as power functions of per-period mean AUC_PAZO and an exponentially declining resistance on the cytotoxic effect. The empirical exponent on capacity growth (n) is fixed at 0.5 for the clinical fit (vs 1 in the paired mouse model) to better describe the tumour-regrowth and long-term-antiangiogenic phases observed in patients."
  reference <- paste(
    "Ouerdani A, Struemper H, Suttle AB, Ouellet D, Ribba B.",
    "Preclinical modeling of tumor growth and angiogenesis inhibition",
    "to describe pazopanib clinical effects in renal cell carcinoma.",
    "CPT Pharmacometrics Syst Pharmacol. 2015;4(11):660-668.",
    "doi:10.1002/psp4.12001.",
    "Erratum/revised version published online 2015-11-12 (Table 1 was replaced;",
    "the on-disk PDF is the corrected version).",
    sep = " "
  )
  vignette <- "Ouerdani_2015_pazopanib"
  units <- list(
    time          = "day",
    dosing        = "n/a (no drug-dosing events; pazopanib exposure enters as the per-period AUC_PAZO covariate, not via a PK ODE)",
    concentration = "mm (sum of longest diameters of target lesions per RECIST 1.1; not a drug concentration)"
  )

  covariateData <- list(
    AUC_PAZO = list(
      description        = "Per-period mean AUC of pazopanib driving the antiangiogenic and cytotoxic drug-effect rates in the Ouerdani 2015 clinical TGI model.",
      units              = "ug*h/mL (= mg*h/L)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying step-wise. Ouerdani 2015 derives mean AUC at each patient's dose level from an Emax fit to pooled mean AUC values reported in five prior pazopanib trials across the 5 mg to 2000 mg daily-dose range (Methods, clinical-data section). The reported population mean was 771.6 ug*h/mL (range 629.4-802.4) corresponding to a mean dose of 727 mg/day (range 473-800) across the 47 RCC patients. Set to 0 in off-treatment periods; the model() block gates both drug-effect rates on AUC_PAZO > 0 so that off-treatment intervals follow the pure-growth dynamics.",
      source_name        = "AUC"
    ),
    TUM_SLD = list(
      description        = "Per-patient observed baseline RECIST 1.1 sum of longest diameters of target lesions; used as the per-subject initial condition for the tumor_size ODE state.",
      units              = "mm",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Ouerdani 2015 sets P0 to the observed baseline SLD per patient so only six structural parameters are estimated (paper Methods, clinical-data section: 'P0 was set to the observed value'). RECIST 1.1 sum of longest diameters of target lesions, in millimetres.",
      source_name        = "SLD0"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 47L,
    n_studies      = 1L,
    age_range      = "43-79 years",
    weight_range   = "not reported in the modelling paper",
    sex_female_pct = NA_real_,
    race_ethnicity = NA,
    disease_state  = "advanced and/or metastatic renal-cell carcinoma of predominantly clear-cell histology; measurable disease by RECIST",
    dose_range     = "800 mg pazopanib once daily, reduced in case of intolerance; mean dose 727 mg/day (range 473-800)",
    regions        = "multicenter open-label Phase 2 study NCT00244764 (regions not detailed in the modelling paper)",
    study          = "Subset of 47 patients from the multicenter, open-label Phase 2 trial NCT00244764 (referenced via Ouerdani 2015 citation 16); the modelling paper used this 47-patient subset.",
    notes          = "Eligibility per Ouerdani 2015 Methods clinical-data section: ECOG PS 0 or 1, adequate haematologic / hepatic / renal function, treatment-naive or single prior systemic immunotherapy with cytokines and/or prior surgery (nephrectomy) and/or radiotherapy. Disease assessments via CT or MRI scheduled at baseline, weeks 8 and 12, then every 8 weeks until progression by RECIST 1.1. No dose interruptions were reported within the analyzed dataset despite the dose-reduction protocol."
  )

  ini({
    # Structural parameters -- clinical column of Ouerdani 2015 Table 1
    # (corrected November 2015 version on disk).
    lk_tumor  <- log(0.0021);  label("Tumour growth rate k (1/day)")                              # Ouerdani 2015 Table 1 clinical k = 0.0021 (RSE 6%)
    lk_cap    <- log(0.0392);  label("Carrying-capacity rate constant b (1/day)")                  # Ouerdani 2015 Table 1 clinical b = 0.0392 (RSE 22%); IIV fixed to 0 per Table 1
    lk_aa0    <- log(0.0023);  label("Baseline antiangiogenic effect rate c0 (1/day)")             # Ouerdani 2015 Table 1 clinical c = 0.0023 (RSE 9%)
    lk_cyto0  <- log(0.0032);  label("Baseline cytotoxic effect rate a0 (1/day)")                  # Ouerdani 2015 Table 1 clinical a = 0.0032 (RSE 2%)
    lk_res    <- log(0.0153);  label("Cytotoxic resistance rate d (1/day; multiplies exp(-d*t))")  # Ouerdani 2015 Table 1 clinical d = 0.0153 (RSE 3%)
    lK0       <- log(329);     label("Initial carrying capacity K0 (mm; SLD scale)")               # Ouerdani 2015 Table 1 clinical K0 = 329 (RSE 25%); IIV fixed to 0

    # Drug-exposure covariate effects on the two drug-effect rates (both
    # estimated in the clinical fit, vs only b_c estimated in the mouse fit).
    e_auc_pazo_k_aa   <- 0.142;   label("Power exponent of AUC_PAZO on antiangiogenic rate c")  # Ouerdani 2015 Table 1 clinical b_c = 0.142 (RSE 7%)
    e_auc_pazo_k_cyto <- 0.125;   label("Power exponent of AUC_PAZO on cytotoxic rate a")        # Ouerdani 2015 Table 1 clinical b_a = 0.125 (RSE 14%)

    # Empirical structural exponent on tumour-driven capacity growth
    # dK/dt = b * P^n - c * K. Fixed at 0.5 for the clinical fit (vs 1 for
    # the mouse fit); paper Results: "For this analysis, it was necessary
    # to apply a correction on the empirical parameter n, which was set in
    # this case to the value 0.5. This value allows the model to better
    # predict the tumor regrowth, as well as the second decrease due to
    # the long-term antiangiogenic effect".
    n <- fixed(0.5); label("Empirical exponent on tumour volume in capacity ODE")                # Ouerdani 2015 Results: clinical fit uses n = 0.5

    # Inter-individual variability (IIV) -- Ouerdani 2015 Table 1 footnote
    # reports IIV as the omega (square root of NONMEM variance) expressed
    # as a percentage. Variance entered as omega^2 = (IIV / 100)^2.
    etalk_tumor ~ 0.6724   # (0.82)^2; clinical IIV on k = 82% (RSE 35%)
    etalk_aa0   ~ 0.0961   # (0.31)^2; clinical IIV on c0 = 31% (RSE 51%)
    etalk_cyto0 ~ 0.3844   # (0.62)^2; clinical IIV on a0 = 62% (RSE 29%)
    etalk_res   ~ 1.0201   # (1.01)^2; clinical IIV on d  = 101% (RSE 45%)
    # No etalK0:    Ouerdani 2015 Table 1 fixes IIV on K0 to 0 (clinical fit).
    # No etalk_cap: Ouerdani 2015 Table 1 fixes IIV on b to 0.

    # Residual error -- combined additive + proportional on SLD. Paper
    # footnote on Table 1: e2 in the clinical column is in mm.
    propSd <- 0.08;  label("Proportional residual SD on SLD (fraction)")                          # Ouerdani 2015 Table 1 clinical e1 = 8% (RSE 2%)
    addSd  <- 1;     label("Additive residual SD on SLD (mm)")                                    # Ouerdani 2015 Table 1 clinical e2 = 1 (RSE 3%)
  })

  model({
    # ----- Individual structural parameters (lognormal IIV where present) -----
    k_tumor <- exp(lk_tumor + etalk_tumor)
    k_cap   <- exp(lk_cap)                                  # no IIV (Table 1)
    a0      <- exp(lk_cyto0 + etalk_cyto0)
    c0      <- exp(lk_aa0   + etalk_aa0)
    k_res   <- exp(lk_res   + etalk_res)
    K0      <- exp(lK0)                                     # no IIV (Table 1)

    # ----- AUC-driven drug-effect rates -----
    # Ouerdani 2015 Methods Equations 2 and 3 (same form as the mouse fit;
    # both exponents are estimated in the clinical fit):
    #     a = a0 * AUC^b_a       (cytotoxic effect rate, 1/day)
    #     c = c0 * AUC^b_c       (antiangiogenic effect rate, 1/day)
    # Both drug-effect rates are gated on AUC_PAZO > 0 so that
    # off-treatment intervals (e.g. permanent dose interruption) carry
    # zero drug effect rather than the AUC^0 = 1 / NONMEM-convention
    # spurious activation.
    if (AUC_PAZO > 0) {
      cyto_rate <- a0 * AUC_PAZO^e_auc_pazo_k_cyto
      aa_rate   <- c0 * AUC_PAZO^e_auc_pazo_k_aa
    } else {
      cyto_rate <- 0
      aa_rate   <- 0
    }

    # ----- ODE system (Ouerdani 2015 Equation 1, clinical fit with n = 0.5) -----
    #   dP/dt = k * P * (1 - P/K) - a * exp(-d*t) * P
    #   dK/dt = b * P^n           - c * K
    # P = tumor_size         (mm, SLD).
    # K = carrying_capacity  (mm, SLD scale); biologically the maximum
    #     tumour size sustainable by the current vasculature.
    # t is rxode2 simulation time in days; treatment start is taken as
    # t = 0 so the exp(-d*t) resistance term decays from 1 at study start.
    d/dt(tumor_size)        <- k_tumor * tumor_size * (1 - tumor_size / carrying_capacity) -
                              cyto_rate * exp(-k_res * t) * tumor_size
    d/dt(carrying_capacity) <- k_cap * tumor_size^n - aa_rate * carrying_capacity

    # Initial state. P0 is per-patient (covariate TUM_SLD); K0 is the
    # population-typical-value parameter (no IIV in the clinical fit).
    tumor_size(0)        <- TUM_SLD
    carrying_capacity(0) <- K0

    # ----- Observation and combined residual error -----
    tumor_size ~ add(addSd) + prop(propSd)
  })
}
