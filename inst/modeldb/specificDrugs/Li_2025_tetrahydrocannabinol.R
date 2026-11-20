Li_2025_tetrahydrocannabinol <- function() {
  description <- paste(
    "Semi-mechanistic parent-metabolite population PK model for oral and",
    "intravenous delta-9-tetrahydrocannabinol (THC) and its sequential",
    "metabolites 11-hydroxy-THC (11-OH-THC) and 11-nor-9-carboxy-THC",
    "(THC-COOH) in adult occasional cannabis users, fitted to digitised mean",
    "concentration-time profiles pooled from 10 published studies. THC and",
    "11-OH-THC each have a central and two peripheral compartments, the",
    "second of which exchanges with an adipose-tissue compartment through a",
    "permeability-surface-area product and an adipose:plasma partition",
    "coefficient; THC-COOH is two-compartment. Oral THC is absorbed first",
    "order from a depot with bioavailability F, all THC clearance forms",
    "11-OH-THC, and a fraction Fm of 11-OH-THC clearance forms THC-COOH",
    "(the remainder excreted unchanged). Variability is between-STUDY, not",
    "between-subject, because the data were study means."
  )
  reference <- paste(
    "Li P, An G. Evaluation of Cannabis Per Se Laws: A Semi-Mechanistic",
    "Pharmacometrics Model for Quantitative Characterization of THC and",
    "Metabolites in Oral Users. J Clin Pharmacol. 2025;65(5):535.",
    "doi:10.1002/jcph.6181.",
    sep = " "
  )
  vignette <- "Li_2025_tetrahydrocannabinol"

  # Methods 'Data Source': 'due to the minimal differences in molecular
  # weight observed between THC and its metabolites (less than 5% per
  # metabolism step), the final dataset opted to use weight-based
  # concentrations rather than molarity.' Formation of each metabolite is
  # therefore 1:1 by mass, dose is in mg and volumes in L, so amount /
  # volume is mg/L and the * 1000 in model() reports ng/mL.
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")
  dosing <- c("depot", "central")

  # No covariate was tested: the data were study-level mean profiles, so
  # no individual covariates were available.
  covariateData <- list()

  compartmentData <- list(
    depot = list(
      analyte = "delta-9-tetrahydrocannabinol",
      units = "mg",
      specimen = "not applicable",
      verified = TRUE
    ),
    central = list(
      analyte = "delta-9-tetrahydrocannabinol",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "delta-9-tetrahydrocannabinol",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    peripheral2 = list(
      analyte = "delta-9-tetrahydrocannabinol",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    adipose = list(
      analyte = "delta-9-tetrahydrocannabinol",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    central_11oh = list(
      analyte = "11-hydroxy-delta-9-tetrahydrocannabinol",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1_11oh = list(
      analyte = "11-hydroxy-delta-9-tetrahydrocannabinol",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    peripheral2_11oh = list(
      analyte = "11-hydroxy-delta-9-tetrahydrocannabinol",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    adipose_11oh = list(
      analyte = "11-hydroxy-delta-9-tetrahydrocannabinol",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    central_thccooh = list(
      analyte = "11-nor-9-carboxy-delta-9-tetrahydrocannabinol",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1_thccooh = list(
      analyte = "11-nor-9-carboxy-delta-9-tetrahydrocannabinol",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = 10L,
    n_dose_groups = 25L,
    data_type = paste(
      "Digitised mean concentration-time profiles (Engauge Digitizer 12.1),",
      "one profile per study arm and analyte; individual data were not",
      "available, so inter-individual variability could not be estimated and",
      "the random effects are inter-STUDY."
    ),
    disease_state = paste(
      "Adults receiving controlled THC doses; occasional cannabis users only",
      "(the frequent-user arms of Newmeyer 2016 and 2017 were excluded) and",
      "fasted conditions only (the fed arm of Lunn 2019 was excluded)."
    ),
    dose_range = paste(
      "Single oral THC 5-90 mg (dronabinol capsules, cannabis decoction,",
      "brownies, commercial edibles, hemp oil) in seven studies and single",
      "intravenous THC 1.6-5 mg in three studies."
    ),
    sampling = "Sampling spanned 4 to 54 h after dose depending on the study.",
    matrices = paste(
      "Plasma in six studies; whole blood in four (Menetrey 2005, Newmeyer",
      "2016, Newmeyer 2017, Spindle 2020), converted to plasma by the authors",
      "with blood:plasma ratios of 0.71 (THC), 0.73 (11-OH-THC) and 0.65",
      "(THC-COOH)."
    ),
    regions = "Pooled published studies (USA, Sweden, Switzerland).",
    notes = paste(
      "Supplement Table S2 lists the 25 study arms (Ohlsson 1980, Menetrey",
      "2005, Lile 2013, Newmeyer 2016, Newmeyer 2017, Lunn 2019, Spindle 2020,",
      "Ewell 2021, Naef 2004, Meyer 2018). Fitted in NONMEM 7.4.3 with FOCEI",
      "and ADVAN14. External validation against Goodwin 2006 (0.13-4.93 mg",
      "THC TID for 5 days; Figure 4)."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # THC. All values are Table 1 'Estimated Parameters from the Final
    # Population PK Model'. Table 1 prints Fm, F and every disposition
    # parameter on the linear scale with no uncertainty interval beyond
    # an RSE%; none is flagged as fixed.
    # ------------------------------------------------------------------
    lfdepot <- log(0.144)
    label("Oral bioavailability of THC, F (fraction)") # Table 1 F = 0.144 (RSE 1%)

    lka <- log(1.46)
    label("First-order absorption rate constant, Ka (1/h)") # Table 1 Ka = 1.46 h-1 (RSE 0.1%)

    lcl <- log(87)
    label("THC clearance, all forming 11-OH-THC, CL_THC (L/h)") # Table 1 CL_THC = 87 L/h (RSE 2.6%)

    lvc <- log(204)
    label("THC central volume, V_THC (L)") # Table 1 V_THC = 204 L (RSE 186.3%)

    lq <- log(59.9)
    label("THC intercompartmental flow to peripheral 1, Q1 (L/h)") # Table 1 Q1 = 59.9 L/h (RSE 0.9%)

    lvp <- log(30.2)
    label("THC peripheral 1 volume, V_THC,peri1 (L)") # Table 1 V_THC,peri1 = 30.2 L (RSE 0.3%)

    lq2 <- log(29.7)
    label("THC intercompartmental flow to peripheral 2, Q2 (L/h)") # Table 1 Q2 = 29.7 L/h (RSE <0.1%)

    lvp2 <- log(88.3)
    label("THC peripheral 2 volume, V_THC,peri2 (L)") # Table 1 V_THC,peri2 = 88.3 L (RSE 8.2%)

    lps_adipose <- log(5.5)
    label("THC permeability-surface-area product, peripheral 2 to adipose, PS (L/h)") # Table 1 PS = 5.5 L/h (RSE <0.1%)

    lv_adipose <- log(91.8)
    label("THC adipose-pool volume, V_THC,adipose (L)") # Table 1 V_THC,adipose = 91.8 L (RSE <0.1%)

    lkp_adipose <- log(14)
    label("THC adipose:plasma partition coefficient, KP (unitless)") # Table 1 KP = 14 (RSE 0.1%)

    # ------------------------------------------------------------------
    # 11-OH-THC (Table 1, second block).
    # ------------------------------------------------------------------
    lcl_11oh <- log(80.4)
    label("11-OH-THC clearance, CL_11-OH-THC (L/h)") # Table 1 CL_11-OH-THC = 80.4 L/h (RSE 0.1%)

    lvc_11oh <- log(49.1)
    label("11-OH-THC central volume, V_11-OH-THC (L)") # Table 1 V_11-OH-THC = 49.1 L (RSE 1.2%)

    lq_11oh <- log(19.9)
    label("11-OH-THC intercompartmental flow to peripheral 1, Q3 (L/h)") # Table 1 Q3 = 19.9 L/h (RSE <0.1%)

    lvp_11oh <- log(9.95)
    label("11-OH-THC peripheral 1 volume, V_11-OH-THC,peri1 (L)") # Table 1 V_11-OH-THC,peri1 = 9.95 L (RSE <0.1%)

    lq2_11oh <- log(47.3)
    label("11-OH-THC intercompartmental flow to peripheral 2, Q4 (L/h)") # Table 1 Q4 = 47.3 L/h (RSE 0.7%)

    lvp2_11oh <- log(84.3)
    label("11-OH-THC peripheral 2 volume, V_11-OH-THC,peri2 (L)") # Table 1 V_11-OH-THC,peri2 = 84.3 L (RSE 1.8%)

    lps_adipose_11oh <- log(5.11)
    label("11-OH-THC permeability-surface-area product, peripheral 2 to adipose, PS2 (L/h)") # Table 1 PS2 = 5.11 L/h (RSE <0.1%)

    lv_adipose_11oh <- log(98.7)
    label("11-OH-THC adipose-pool volume, V_11-OH-THC,adipose (L)") # Table 1 V_11-OH-THC,adipose = 98.7 L (RSE <0.1%)

    lkp_adipose_11oh <- log(9.87)
    label("11-OH-THC adipose:plasma partition coefficient, KP2 (unitless)") # Table 1 KP2 = 9.87 (RSE <0.1%)

    fm <- 0.904
    label("Fraction of 11-OH-THC clearance forming THC-COOH, Fm (fraction)") # Table 1 Fm = 0.904 (RSE 0.6%)

    # ------------------------------------------------------------------
    # THC-COOH (Table 1, third block).
    # ------------------------------------------------------------------
    lcl_thccooh <- log(4.34)
    label("THC-COOH clearance, CL_THC-COOH (L/h)") # Table 1 CL_THC-COOH = 4.34 L/h (RSE 0.3%)

    lvc_thccooh <- log(0.986)
    label("THC-COOH central volume, V_THC-COOH (L)") # Table 1 V_THC-COOH = 0.986 L (RSE 0.1%)

    lq_thccooh <- log(6.35)
    label("THC-COOH intercompartmental flow, Q5 (L/h)") # Table 1 Q5 = 6.35 L/h (RSE <0.1%)

    lvp_thccooh <- log(27.6)
    label("THC-COOH peripheral volume, V_THC-COOH,peri (L)") # Table 1 V_THC-COOH,peri = 27.6 L (RSE 0.1%)

    # ------------------------------------------------------------------
    # Inter-STUDY variability, exponential (Equation 12). Table 1 reports
    # the variances omega^2 directly. Results 'Parameter Estimation'
    # names the four parameters as 'CL, V2, V3, and V4'.
    # ------------------------------------------------------------------
    etalcl ~ 0.506 # Table 1 'omega2 CL_THC' = 0.506 (shrinkage 7.2%)
    etalvc ~ 1.19 # Table 1 'omega2 V_THC' = 1.19 (shrinkage 0.1%)
    etalvc_11oh ~ 1.08 # Table 1 'omega2 V_11-OH-THC' = 1.08 (shrinkage 0.1%)
    etalvc_thccooh ~ 0.846 # Table 1 'omega2 V_THC-COOH' = 0.846 (shrinkage 37.5%)

    # ------------------------------------------------------------------
    # Residual error. Table 1 reports only a proportional VARIANCE per
    # analyte, so propSd = sqrt(sigma^2). Methods Equation 13 describes a
    # combined proportional + additive model as the one evaluated, but no
    # additive component is reported for the final model.
    # ------------------------------------------------------------------
    propSd <- sqrt(0.375)
    label("Proportional residual error for THC (fraction)") # Table 1 'sigma2 prop,THC' = 0.375; SD = sqrt(0.375) = 0.612

    propSd_11oh <- sqrt(0.31)
    label("Proportional residual error for 11-OH-THC (fraction)") # Table 1 'sigma2 prop,11-OH-THC' = 0.31; SD = 0.557

    propSd_thccooh <- sqrt(0.144)
    label("Proportional residual error for THC-COOH (fraction)") # Table 1 'sigma2 prop,THC-COOH' = 0.144; SD = 0.379
  })
  model({
    # ---- Individual parameters ---------------------------------------
    ka <- exp(lka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp)
    q2 <- exp(lq2)
    vp2 <- exp(lvp2)
    ps_adipose <- exp(lps_adipose)
    v_adipose <- exp(lv_adipose)
    kp_adipose <- exp(lkp_adipose)

    cl_11oh <- exp(lcl_11oh)
    vc_11oh <- exp(lvc_11oh + etalvc_11oh)
    q_11oh <- exp(lq_11oh)
    vp_11oh <- exp(lvp_11oh)
    q2_11oh <- exp(lq2_11oh)
    vp2_11oh <- exp(lvp2_11oh)
    ps_adipose_11oh <- exp(lps_adipose_11oh)
    v_adipose_11oh <- exp(lv_adipose_11oh)
    kp_adipose_11oh <- exp(lkp_adipose_11oh)

    cl_thccooh <- exp(lcl_thccooh)
    vc_thccooh <- exp(lvc_thccooh + etalvc_thccooh)
    q_thccooh <- exp(lq_thccooh)
    vp_thccooh <- exp(lvp_thccooh)

    # ---- Concentrations (mg/L) ---------------------------------------
    c_central <- central / vc
    c_peri1 <- peripheral1 / vp
    c_peri2 <- peripheral2 / vp2
    c_adipose <- adipose / v_adipose
    c_central_11oh <- central_11oh / vc_11oh
    c_peri1_11oh <- peripheral1_11oh / vp_11oh
    c_peri2_11oh <- peripheral2_11oh / vp2_11oh
    c_adipose_11oh <- adipose_11oh / v_adipose_11oh
    c_central_thccooh <- central_thccooh / vc_thccooh
    c_peri1_thccooh <- peripheral1_thccooh / vp_thccooh

    # ---- ODEs: Equations 1-11 ----------------------------------------
    # Adipose exchange (Eqs. 6-7, 9-10) is permeability-limited between
    # the second peripheral compartment and adipose, with the adipose
    # concentration divided by the partition coefficient.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - (cl + q + q2) * c_central +
      q * c_peri1 + q2 * c_peri2
    d/dt(peripheral1) <- q * (c_central - c_peri1)
    d/dt(peripheral2) <- q2 * (c_central - c_peri2) -
      ps_adipose * (c_peri2 - c_adipose / kp_adipose)
    d/dt(adipose) <- ps_adipose * (c_peri2 - c_adipose / kp_adipose)

    d/dt(central_11oh) <- cl * c_central -
      (cl_11oh + q_11oh + q2_11oh) * c_central_11oh +
      q_11oh * c_peri1_11oh + q2_11oh * c_peri2_11oh
    d/dt(peripheral1_11oh) <- q_11oh * (c_central_11oh - c_peri1_11oh)
    d/dt(peripheral2_11oh) <- q2_11oh * (c_central_11oh - c_peri2_11oh) -
      ps_adipose_11oh * (c_peri2_11oh - c_adipose_11oh / kp_adipose_11oh)
    d/dt(adipose_11oh) <- ps_adipose_11oh *
      (c_peri2_11oh - c_adipose_11oh / kp_adipose_11oh)

    # (1 - fm) of 11-OH-THC clearance is renal excretion of unchanged
    # 11-OH-THC and leaves the system.
    d/dt(central_thccooh) <- cl_11oh * fm * c_central_11oh -
      (cl_thccooh + q_thccooh) * c_central_thccooh +
      q_thccooh * c_peri1_thccooh
    d/dt(peripheral1_thccooh) <- q_thccooh *
      (c_central_thccooh - c_peri1_thccooh)

    # Bioavailability applies to the oral depot only; intravenous doses
    # go to central with F = 1.
    f(depot) <- exp(lfdepot)

    # ---- Observations (ng/mL) ----------------------------------------
    Cc <- 1000 * c_central
    Cc_11oh <- 1000 * c_central_11oh
    Cc_thccooh <- 1000 * c_central_thccooh

    Cc ~ prop(propSd)
    Cc_11oh ~ prop(propSd_11oh)
    Cc_thccooh ~ prop(propSd_thccooh)
  })
}
