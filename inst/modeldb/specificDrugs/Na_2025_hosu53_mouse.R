Na_2025_hosu53_mouse <- function() {
  description <- paste(
    "Preclinical (mouse). Coupled two-compartment population PK and",
    "indirect-response turnover PD model for the oral dihydroorotate",
    "dehydrogenase (DHODH) inhibitor HOSU-53 (JBZ-001) and its plasma",
    "pharmacodynamic biomarker dihydroorotate (DHO) in mice",
    "(Na 2025 Supplementary Table S3). PK is two-compartment with",
    "first-order absorption, logit-normal bioavailability, and linear",
    "elimination. DHO is described by an indirect-response turnover model",
    "in which HOSU-53 inhibits the first-order degradation of DHO through a",
    "sigmoidal Imax function with Imax fixed at 1:",
    "dDHO/dt = kin - kout * (1 - Cc^hill / (Cc^hill + ic50^hill)) * DHO,",
    "with kin = rbase * kout and DHO(0) = rbase (Na 2025 Section 3.4).",
    "PK and PD parameters were estimated simultaneously for the mouse data",
    "(Na 2025 Section 2.6). Dose level, sex, salt form, and vehicle",
    "formulation were screened as covariates and none were retained."
  )
  reference <- paste(
    "Na JY, Hai M, Kim K, Vibhute SM, Bennett CE, Coss CC, Phelps MA.",
    "Translational Pharmacokinetic-Pharmacodynamic Modeling of a Novel Oral",
    "Dihydroorotate Dehydrogenase (DHODH) Inhibitor, HOSU-53 (JBZ-001).",
    "Pharmaceutics. 2025;17(4):412. doi:10.3390/pharmaceutics17040412"
  )
  vignette <- "Na_2025_hosu53"
  units <- list(time = "h", dosing = "nmol", concentration = "umol/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Amount units are nmol because Na 2025 reports the
  # volumes in mL (Supplementary Table S3) and both concentrations
  # (HOSU-53, DHO) in umol/L: nmol / mL = umol/L. The paper does not state
  # the dose-amount unit used in the Monolix dataset and does not report the
  # HOSU-53 molecular weight, so mg/kg doses cannot be converted here; see
  # the vignette "Assumptions and deviations" section.
  compartmentData <- list(
    depot       = list(analyte = "HOSU-53", units = "nmol", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "HOSU-53", units = "nmol", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "HOSU-53", units = "nmol", specimen = "plasma", verified = TRUE),
    dho         = list(analyte = "dihydroorotate", units = "umol/L", specimen = "plasma", verified = TRUE)
  )

  # Na 2025 Section 2.7 screened dose level, sex, salt form, and vehicle
  # formulation on the PK and PD parameters by stepwise selection
  # (forward dOFV 3.84, backward dOFV 5.99). Section 3.3: "none of the
  # covariates tested on the PK parameters were found to be significant".
  # They are documented here rather than in covariateData because the final
  # model references none of them.
  covariatesDataExcluded <- list(
    SEXF = list(
      description        = "Female sex indicator. Screened on the PK and PD parameters (Na 2025 Section 2.7); not retained.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Both sexes contributed to the mouse dataset (Na 2025 Section 2.2: male mice at 4, 10, 20, and 30 mg/kg; female mice at 10 mg/kg). Section 3.3 and the Discussion report that sex was not a significant covariate."
    ),
    DOSE_HOSU53_MGKG = list(
      description        = "Administered HOSU-53 dose level in mg/kg. Screened as a categorical covariate on the PK and PD parameters (Na 2025 Section 2.7); not retained.",
      units              = "mg/kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Dose levels contributing to the modelled mouse dataset were 3 mg/kg (IV), 4 mg/kg, and 10 mg/kg (Na 2025 Section 3.1). Section 3.2 reports less-than-dose-proportional exposure in mice over the wider 4-30 mg/kg range, which the linear-elimination structure retained here does not reproduce."
    ),
    FORM_HOSU53_SALT = list(
      description        = "HOSU-53 salt form administered (sodium salt vs lysine salt). Screened on the PK parameters (Na 2025 Section 2.7); not retained.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (sodium salt)",
      notes              = "Na 2025 Section 2.1 reports bridging PK studies between the sodium and lysine salts that 'demonstrated near-identical results'; the covariate screen confirmed no significant effect."
    ),
    FORM_HOSU53_VEHICLE = list(
      description        = "Dosing vehicle used for the oral formulation (40% or 20% HPBCD). Screened on the PK parameters (Na 2025 Section 2.7); not retained.",
      units              = "(categorical)",
      type               = "categorical",
      reference_category = "40% hydroxypropyl-beta-cyclodextrin (HPBCD)",
      notes              = "The second mouse single-dose study used 20% HPBCD specifically 'to reduce the potential impact of formulation components on HOSU-53 disposition' (Na 2025 Section 2.2); the covariate screen found no significant vehicle effect."
    )
  )

  population <- list(
    species          = "mouse (non-tumour-bearing, and male and female triple-immunodeficient NCG mice bearing MOLM-13 disseminated AML xenografts)",
    n_subjects       = 45L,
    n_studies        = 3L,
    n_pk_samples     = 309L,
    n_pd_samples     = 253L,
    sex              = "male and female",
    disease_state    = "healthy mice (single-dose PK/PD studies) and MOLM-13 human disseminated acute myeloid leukaemia xenograft-bearing NCG mice (repeat-dose study)",
    dose_range       = "3 mg/kg IV and 4 and 10 mg/kg PO in the modelled dataset (3/4/10 mg/kg = 9/9/27 mice). The broader mouse programme also dosed 20 mg/kg PO QD and 30 mg/kg PO twice weekly, which did not contribute to the modelled dataset. All PO dosing was by gavage and all IV dosing was by bolus.",
    blq_handling     = "No HOSU-53 PK samples were below the limit of quantification (excluding pre-dose). 9 of 253 DHO PD samples (<4%) were below the limit of quantification and were excluded from the analysis (Na 2025 Section 3.1).",
    assay_range      = "UHPLC-MS/MS calibration range 10.0-5000 ng/mL for HOSU-53 and 10.0-30,000 ng/mL for DHO in mouse plasma (Na 2025 Supplementary Materials).",
    notes            = "Study-level detail is in Na 2025 Supplementary Table S1. Modelling was performed with the SAEM algorithm in Monolix 2024R1; PK and PD parameters were estimated simultaneously (Section 2.6) and parameter precision was assessed by nonparametric bootstrap. Section 3.2 notes that mouse elimination was slower than in rats and dogs (mean NCA half-life 30 h vs 6-8 h), possibly because of higher plasma protein binding."
  )

  ini({
    # ------------------------------------------------------------------
    # PK structural parameters (Na 2025 Supplementary Table S3).
    #
    # Volumes and clearances are reported in mL and mL/h. Concentrations
    # are umol/L (Table S3 units for R0 and IC50; Supplementary Figure S1
    # axis labels), so the model amount unit is nmol (nmol / mL = umol/L).
    #
    # Table S3 labels the clearance row "CL/F", but F is estimated
    # separately in the same table and the mouse fit used pooled IV
    # (3 mg/kg) and PO data, so the value is a true (not apparent)
    # clearance. Dimensional check: the resulting terminal half-life is
    # log(2) / beta = 24 h, consistent with the 30 h mean NCA half-life
    # reported for mice in Section 3.2.
    # ------------------------------------------------------------------
    logitfdepot <- log(0.53 / (1 - 0.53)); label("Logit of oral bioavailability F (unitless; F = 0.53)")  # Table S3 (F = 0.53); logit scale per Section 2.5
    lka         <- log(1.0);               label("First-order absorption rate constant (Ka, 1/h)")        # Table S3 (Ka = 1.0 /h)
    lcl         <- log(0.07);              label("Clearance (CL, mL/h)")                                  # Table S3 (CL = 0.07 mL/h, RSE 5.0%)
    lvc         <- log(1.2);               label("Central volume of distribution (V1, mL)")               # Table S3 (V1 = 1.2 mL, RSE 8.7%)
    # Table S3's row header reads "Q (mL/h)" but its abbreviation footnote reads
    # "Q, intercompartment clearance (L/h)". mL/h is the correct reading: it is
    # the unit consistent with CL (0.07 mL/h) and V1 (1.2 mL) in the same table,
    # and only mL/h reproduces the ~30 h NCA half-life of Section 3.2.
    lq          <- log(1.6);               label("Intercompartmental clearance (Q, mL/h)")                # Table S3 (Q = 1.6 mL/h, RSE 1.7%)
    lvp         <- log(1.2);               label("Peripheral volume of distribution (V2, mL)")            # Table S3 (V2 = 1.2 mL, RSE 5.3%)

    # ------------------------------------------------------------------
    # DHO indirect-response (turnover) PD parameters (Na 2025 Table S3).
    #
    # Section 3.4: dR/dt = kin - kout * (1 - Cp^gamma/(Cp^gamma + IC50^gamma)) * R
    # with kin = R0 * kout and R0 the baseline (steady-state) plasma DHO.
    # The printed final-model equation carries no Imax term, i.e. Imax = 1;
    # the general form with Imax appears in Section 2.11 (human prediction).
    # The paper's "gamma (sigmoidicity of the drug effect)" is the Hill
    # exponent of the sigmoidal Imax function, so it is the canonical
    # `lhill` rather than `lgamma`.
    # ------------------------------------------------------------------
    lrbase <- log(0.01); label("Baseline (steady-state) plasma DHO concentration (R0, umol/L)")            # Table S3 (R0 = 0.01 umol/L, RSE 1.4%)
    lkout  <- log(155);  label("First-order DHO degradation rate constant (Kout, 1/h)")                    # Table S3 (Kout = 155 /h, RSE 12.7%)
    lic50  <- log(1.55); label("HOSU-53 concentration giving half-maximal inhibition of DHO degradation (IC50, umol/L)")  # Table S3 (IC50 = 1.55 umol/L, RSE 8.9%)
    lhill  <- log(1.71); label("Sigmoidicity (Hill) exponent of the inhibitory function (gamma, unitless)")  # Table S3 (gamma = 1.71, RSE 0.1%)
    limax  <- fixed(log(1)); label("Maximum fractional inhibition of DHO degradation (Imax, unitless)")      # Section 3.4 equation carries no Imax term, i.e. Imax = 1

    # ------------------------------------------------------------------
    # Inter-individual variability (Na 2025 Supplementary Table S3).
    #
    # Monolix reports omega on the STANDARD-DEVIATION scale of the random
    # effect (the same convention as the "proportional residual error" rows
    # of this table, which are the SD-scale b parameter of the Monolix
    # error model). nlmixr2 ini() takes the VARIANCE, so each tabulated
    # omega is squared here. Table S3 reports no IIV on F, Ka, or Q.
    # ------------------------------------------------------------------
    etalcl   ~ 0.32^2   # Table S3 (IIV CL = 0.32, RSE 11.6%)
    etalvc   ~ 0.53^2   # Table S3 (IIV V1 = 0.53, RSE 12.1%)
    etalvp   ~ 0.16^2   # Table S3 (IIV V2 = 0.16, RSE 41.8%)
    etalkout ~ 0.62^2   # Table S3 (IIV Kout = 0.62, RSE 12.7%)
    etalic50 ~ 0.53^2   # Table S3 (IIV IC50 = 0.53, RSE 8.9%)

    # ------------------------------------------------------------------
    # Residual error (Na 2025 Table S3). Section 3.3: "Proportional and
    # combined error models were selected for the mouse and dog data,
    # respectively", so the mouse PK error is proportional only.
    # ------------------------------------------------------------------
    propSd     <- 0.18; label("Proportional residual SD for HOSU-53 (fraction)")      # Table S3 (proportional residual error (PK) 0.18, RSE 5.5%)
    propSd_dho <- 0.46; label("Proportional residual SD for plasma DHO (fraction)")   # Table S3 (proportional residual error (PD) 0.46, RSE 6.3%)
  })

  model({
    # --- Individual PK parameters -------------------------------------
    # F is logit-normal (Section 2.5) but carries no IIV in the mouse fit;
    # the remaining parameters are log-normal.
    fdepot <- 1 / (1 + exp(-logitfdepot))
    ka     <- exp(lka)
    cl     <- exp(lcl + etalcl)
    vc     <- exp(lvc + etalvc)
    q      <- exp(lq)
    vp     <- exp(lvp + etalvp)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # --- HOSU-53 disposition ------------------------------------------
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                k12 * central - k21 * peripheral1

    f(depot) <- fdepot

    Cc <- central / vc   # HOSU-53 plasma concentration, nmol/mL = umol/L

    # --- DHO indirect response (degradation inhibition) ---------------
    # Na 2025 Section 3.4:
    #   dR/dt = kin - kout * (1 - Cp^gamma / (Cp^gamma + IC50^gamma)) * R
    #   kin   = R0 * kout,  R(0) = R0
    rbase <- exp(lrbase)
    kout  <- exp(lkout + etalkout)
    ic50  <- exp(lic50 + etalic50)
    hill  <- exp(lhill)
    imax  <- exp(limax)

    kin_dho <- rbase * kout
    inh_dho <- imax * Cc^hill / (Cc^hill + ic50^hill)

    dho(0)    <- rbase
    d/dt(dho) <- kin_dho - kout * (1 - inh_dho) * dho

    Cc  ~ prop(propSd)
    dho ~ prop(propSd_dho)
  })
}
