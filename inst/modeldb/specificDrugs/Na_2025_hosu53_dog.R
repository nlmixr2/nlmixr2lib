Na_2025_hosu53_dog <- function() {
  description <- paste(
    "Preclinical (beagle dog). Coupled two-compartment population PK and",
    "indirect-response turnover PD model for the oral dihydroorotate",
    "dehydrogenase (DHODH) inhibitor HOSU-53 (JBZ-001) and its plasma",
    "pharmacodynamic biomarker dihydroorotate (DHO) in beagle dogs",
    "(Na 2025 Table 3). PK is two-compartment with first-order absorption,",
    "logit-normal bioavailability, and linear elimination. DHO is described",
    "by an indirect-response turnover model in which HOSU-53 inhibits the",
    "first-order degradation of DHO through a sigmoidal Imax function with",
    "Imax fixed at 1:",
    "dDHO/dt = kin - kout * (1 - Cc^hill / (Cc^hill + ic50^hill)) * DHO,",
    "with kin = rbase * kout and DHO(0) = rbase (Na 2025 Section 3.4).",
    "The dog model was fit sequentially: the four disposition parameters",
    "(CL, V1, Q, V2), the IIV on CL, and both PK residual-error terms were",
    "fixed at the IV-only fit so that bioavailability (F) and the absorption",
    "rate constant (Ka) could be estimated from the pooled IV + PO data",
    "(Na 2025 Section 3.3). Dose level, sex, salt form, and vehicle",
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
  # volumes in mL (Table 3) and both concentrations (HOSU-53, DHO) in
  # umol/L: nmol / mL = umol/L. The paper does not state the dose-amount
  # unit used in the Monolix dataset and does not report the HOSU-53
  # molecular weight, so mg/kg doses cannot be converted here; see the
  # vignette "Assumptions and deviations" section.
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
      notes              = "Both sexes contributed to the dog dataset. Section 3.3 and the Discussion report that sex was not a significant covariate on any PK or PD parameter."
    ),
    DOSE_HOSU53_MGKG = list(
      description        = "Administered HOSU-53 dose level in mg/kg. Screened as a categorical covariate on the PK and PD parameters (Na 2025 Section 2.7); not retained.",
      units              = "mg/kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Included dog dose levels were 0.2, 0.3, 0.6, and 3 mg/kg. Na 2025 Section 3.2 reports dose-proportional exposure in dogs over 0.3-3 mg/kg, consistent with the linear-elimination structure retained here."
    ),
    FORM_HOSU53_SALT = list(
      description        = "HOSU-53 salt form administered (sodium salt vs lysine salt). Screened on the PK parameters (Na 2025 Section 2.7); not retained.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (sodium salt)",
      notes              = "Na 2025 Section 2.1 reports bridging PK studies between the sodium and lysine salts that 'demonstrated near-identical results'; the covariate screen confirmed no significant effect."
    ),
    FORM_HOSU53_VEHICLE = list(
      description        = "Dosing vehicle used for the oral formulation (40% HPBCD, deionized water, or 50 mM Tris buffer). Screened on the PK parameters (Na 2025 Section 2.7); not retained.",
      units              = "(categorical)",
      type               = "categorical",
      reference_category = "40% hydroxypropyl-beta-cyclodextrin (HPBCD)",
      notes              = "Vehicle was varied across the dog studies (Na 2025 Section 2.4) and included in the covariate screen because of its expected effect on solubility and absorption; not significant."
    )
  )

  population <- list(
    species          = "beagle dog",
    n_subjects       = 39L,
    n_studies        = 4L,
    n_pk_samples     = 507L,
    n_pd_samples     = 431L,
    sex              = "male and female",
    disease_state    = "healthy beagle dogs (non-GLP PK/PD and GLP toxicokinetic/toxicodynamic studies)",
    dose_range       = "0.2, 0.3, 0.6, and 3 mg/kg PO (0.2/0.3/0.6/3 mg/kg = 8/12/10/9 dogs); an initial single-dose study used 10 mg/kg PO and 3 mg/kg IV. All PO dosing was by gavage and all IV dosing was by bolus.",
    excluded_arms    = "The 1 mg/kg PO QD arm of the 28-day GLP study was excluded from the modelling dataset because of adverse findings (decreased food consumption, lower body weights, intestinal mucosal epithelial cell loss) (Na 2025 Section 3.1).",
    blq_handling     = "No HOSU-53 PK samples were below the limit of quantification (excluding pre-dose). 92 of 431 DHO PD samples (21.3%), including pre-dose samples, were below the limit of quantification and were excluded from the analysis (Na 2025 Section 3.1).",
    assay_range      = "UHPLC-MS/MS calibration range 5.0-5000 ng/mL for HOSU-53 and 10.0-30,000 ng/mL for DHO in dog plasma (Na 2025 Supplementary Materials).",
    noael            = "0.6 mg/kg PO QD in the 28-day GLP repeat-dose study; the highest non-severely toxic dose (HNSTD) was 0.6 mg/kg/day (Na 2025 Sections 3.2 and 4).",
    notes            = "Study-level detail is in Na 2025 Supplementary Table S1. Modelling was performed with the SAEM algorithm in Monolix 2024R1; parameter precision was assessed by nonparametric bootstrap."
  )

  ini({
    # ------------------------------------------------------------------
    # PK structural parameters (Na 2025 Table 3).
    #
    # Volumes and clearances are reported in mL and mL/h. Concentrations
    # are umol/L (Figure 1 axis labels, Table 3 units for R0 and IC50), so
    # the model amount unit is nmol (nmol / mL = umol/L).
    #
    # CL, V1, Q, and V2 were fixed at the IV-only dog fit so that F and Ka
    # could be estimated from the pooled IV + PO dataset (Section 3.3);
    # Table 3 marks each of them "Fix". Table 3 labels the clearance row
    # "CL/F", but F is estimated separately in the same table and the four
    # disposition parameters come from the IV-only fit, so the value is a
    # true (not apparent) clearance. Dimensional check: the resulting
    # terminal half-life is log(2) / beta = 7.1 h, matching the 6-8 h mean
    # NCA half-life reported for dogs in Section 3.2.
    # ------------------------------------------------------------------
    logitfdepot <- log(0.67 / (1 - 0.67)); label("Logit of oral bioavailability F (unitless; F = 0.67)")  # Table 3 (F = 0.67, RSE 2.8%); logit scale per Section 2.5
    lka         <- log(1.53);              label("First-order absorption rate constant (Ka, 1/h)")        # Table 3 (Ka = 1.53 /h, RSE 12.6%)
    lcl         <- fixed(log(150));        label("Clearance (CL, mL/h)")                                  # Table 3 (CL = 150 mL/h, from the IV-only fit)
    lvc         <- fixed(log(980));        label("Central volume of distribution (V1, mL)")               # Table 3 (V1 = 980 mL, from the IV-only fit)
    # Table 3's row header reads "Q (mL/h)" but its abbreviation footnote reads
    # "Q, intercompartment clearance (L/h)". mL/h is the correct reading: it is
    # the unit consistent with CL (150 mL/h) and V1 (980 mL) in the same table,
    # 400 L/h is physiologically impossible in a dog, and only mL/h reproduces
    # the 6-8 h NCA half-life of Section 3.2.
    lq          <- fixed(log(400));        label("Intercompartmental clearance (Q, mL/h)")                # Table 3 (Q = 400 mL/h, from the IV-only fit)
    lvp         <- fixed(log(490));        label("Peripheral volume of distribution (V2, mL)")            # Table 3 (V2 = 490 mL, from the IV-only fit)

    # ------------------------------------------------------------------
    # DHO indirect-response (turnover) PD parameters (Na 2025 Table 3).
    #
    # Section 3.4: dR/dt = kin - kout * (1 - Cp^gamma/(Cp^gamma + IC50^gamma)) * R
    # with kin = R0 * kout and R0 the baseline (steady-state) plasma DHO.
    # The printed final-model equation carries no Imax term, i.e. Imax = 1;
    # the general form with Imax appears in Section 2.11 (human prediction).
    # The paper's "gamma (sigmoidicity of the drug effect)" is the Hill
    # exponent of the sigmoidal Imax function, so it is the canonical
    # `lhill` rather than `lgamma`.
    # ------------------------------------------------------------------
    lrbase <- log(0.06); label("Baseline (steady-state) plasma DHO concentration (R0, umol/L)")            # Table 3 (R0 = 0.06 umol/L, RSE 1.6%)
    lkout  <- log(52);   label("First-order DHO degradation rate constant (Kout, 1/h)")                    # Table 3 (Kout = 52 /h, RSE 11.4%)
    lic50  <- log(0.1);  label("HOSU-53 concentration giving half-maximal inhibition of DHO degradation (IC50, umol/L)")  # Table 3 (IC50 = 0.1 umol/L, RSE 10.3%)
    lhill  <- log(1.9);  label("Sigmoidicity (Hill) exponent of the inhibitory function (gamma, unitless)")  # Table 3 (gamma = 1.9, RSE 1.6%)
    limax  <- fixed(log(1)); label("Maximum fractional inhibition of DHO degradation (Imax, unitless)")      # Section 3.4 equation carries no Imax term, i.e. Imax = 1

    # ------------------------------------------------------------------
    # Inter-individual variability (Na 2025 Table 3).
    #
    # Monolix reports omega on the STANDARD-DEVIATION scale of the random
    # effect (the same convention as the "proportional residual error" rows
    # of this table, which are the SD-scale b parameter of the Monolix
    # combined error model). nlmixr2 ini() takes the VARIANCE, so each
    # tabulated omega is squared here. IIV on F is on the logit scale
    # (Section 2.5: F was modelled with a logit-normal distribution); all
    # other IIVs are log-normal (Section 2.5 exponential model).
    # ------------------------------------------------------------------
    etalogitfdepot ~ 0.24^2         # Table 3 (IIV F = 0.24, RSE 60.1%)
    etalka         ~ 0.59^2         # Table 3 (IIV Ka = 0.59, RSE 20.1%)
    etalcl         ~ fixed(0.2^2)   # Table 3 (IIV CL = 0.2, from the IV-only fit)
    etalkout       ~ 0.42^2         # Table 3 (IIV Kout = 0.42, RSE 14.9%)
    etalic50       ~ 0.44^2         # Table 3 (IIV IC50 = 0.44, RSE 12.9%)

    # ------------------------------------------------------------------
    # Residual error (Na 2025 Table 3). The PK error is the Monolix
    # combined (additive + proportional) model, both terms fixed at the
    # IV-only fit. The PD error is proportional only.
    # ------------------------------------------------------------------
    addSd      <- fixed(0.0033); label("Additive residual SD for HOSU-53 (umol/L)")            # Table 3 (additive residual error 0.0033, from the IV-only fit)
    propSd     <- fixed(0.45);   label("Proportional residual SD for HOSU-53 (fraction)")      # Table 3 (proportional residual error 0.45, from the IV-only fit)
    propSd_dho <- 0.55;          label("Proportional residual SD for plasma DHO (fraction)")   # Table 3 (proportional residual error (PD) 0.55, RSE 5.1%)
  })

  model({
    # --- Individual PK parameters -------------------------------------
    # F is logit-normal (Section 2.5); the remaining parameters are
    # log-normal. CL carries the IIV that was fixed from the IV-only fit.
    logit_f <- logitfdepot + etalogitfdepot
    fdepot  <- 1 / (1 + exp(-logit_f))
    ka     <- exp(lka + etalka)
    cl     <- exp(lcl + etalcl)
    vc     <- exp(lvc)
    q      <- exp(lq)
    vp     <- exp(lvp)

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

    Cc  ~ add(addSd) + prop(propSd)
    dho ~ prop(propSd_dho)
  })
}
