Kim_2022_SCR430_rat_repeatedDosing_pbpk <- function() {
  description <- paste(
    "Preclinical (rat).",
    "PBPK (minimal, IVIVE-based, Berkeley Madonna).",
    "SCR430 (investigational sorafenib derivative) in female Sprague-Dawley",
    "rats pretreated with oral SCR430 30 mg/kg once daily for 14 days",
    "(repeated-dosing group, auto-induced hepatic uptake), after a 2-day",
    "washout and a 3 mg/kg intravenous bolus.",
    "A systemic blood compartment exchanges with a single-adjustment",
    "compartment by first-order rate constants, and the liver is a",
    "five-unit tandem (dispersion-model emulation) of extracellular and",
    "hepatocellular spaces linked in series by hepatic blood flow.",
    "Hepatocyte entry is active uptake plus passive diffusion acting on",
    "the unbound blood concentration; passive efflux and intrinsic",
    "metabolic clearance act on the unbound hepatocellular concentration.",
    "The hepatic uptake, diffusion, metabolism and intracellular-binding",
    "inputs come from hepatocytes isolated from repeated-dosing rats (IVIVE",
    "Method 2) and carry the source paper's empirical 1.6-fold scaling",
    "factor on the permeability clearances; the compartment volume and",
    "exchange rate constants are those fitted for the control group",
    "('Equal to Control').",
    "Between-animal variability on the four hepatic inputs is the source",
    "paper's lognormal Monte Carlo distribution. Companion model:",
    "Kim_2022_SCR430_rat_control_pbpk (untreated rats)."
  )
  reference <- paste(
    "Kim M-C, Lee Y-J.",
    "Analysis of Time-Dependent Pharmacokinetics Using In Vitro-In Vivo",
    "Extrapolation and Physiologically Based Pharmacokinetic Modeling.",
    "Pharmaceutics. 2022;14(12):2562.",
    "doi:10.3390/pharmaceutics14122562.",
    sep = " "
  )
  vignette <- "Kim_2022_SCR430_rat_pbpk"

  # `is_liver<n>` = liver extracellular (sinusoidal) unit n and
  # `int_liver<n>` = hepatocellular unit n of the five-unit tandem liver
  # (source paper Figure 7 'Inlet 1-5' / 'liver 1-5'); the unit index is not
  # part of the registered `is_<organ>` / `int_<organ>` regex, so the ten
  # states are declared here, as in Asaumi_2019_coproporphyrin_I_rifampicin_pbpk.
  paper_specific_compartments <- c(
    "is_liver1",
    "is_liver2",
    "is_liver3",
    "is_liver4",
    "is_liver5",
    "int_liver1",
    "int_liver2",
    "int_liver3",
    "int_liver4",
    "int_liver5"
  )

  units <- list(time = "h", dosing = "mg/kg", concentration = "ug/mL")

  # Every state holds an amount PER KILOGRAM of body weight: the source
  # paper reports every volume, flow and clearance per kg and doses per kg,
  # so body weight cancels out of every concentration.
  compartmentData <- list(
    central = list(analyte = "SCR430", units = "mg/kg", specimen = "whole blood", verified = TRUE),
    peripheral1 = list(analyte = "SCR430", units = "mg/kg", specimen = "whole blood", verified = TRUE),
    is_liver1 = list(analyte = "SCR430", units = "mg/kg", specimen = "whole blood", verified = TRUE),
    is_liver2 = list(analyte = "SCR430", units = "mg/kg", specimen = "whole blood", verified = TRUE),
    is_liver3 = list(analyte = "SCR430", units = "mg/kg", specimen = "whole blood", verified = TRUE),
    is_liver4 = list(analyte = "SCR430", units = "mg/kg", specimen = "whole blood", verified = TRUE),
    is_liver5 = list(analyte = "SCR430", units = "mg/kg", specimen = "whole blood", verified = TRUE),
    int_liver1 = list(analyte = "SCR430", units = "mg/kg", specimen = "tissue", verified = TRUE),
    int_liver2 = list(analyte = "SCR430", units = "mg/kg", specimen = "tissue", verified = TRUE),
    int_liver3 = list(analyte = "SCR430", units = "mg/kg", specimen = "tissue", verified = TRUE),
    int_liver4 = list(analyte = "SCR430", units = "mg/kg", specimen = "tissue", verified = TRUE),
    int_liver5 = list(analyte = "SCR430", units = "mg/kg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "rat (female Sprague-Dawley)",
    n_subjects = 6L,
    n_studies = 1L,
    age_range = "8 weeks",
    weight_range = "190-205 g",
    sex_female_pct = 100,
    disease_state = paste(
      "Healthy rats given SCR430 30 mg/kg in vehicle orally once daily for",
      "14 days, then a 2-day washout before the pharmacokinetic study",
      "(source paper Figure 1, black arm)."
    ),
    dose_range = paste(
      "SCR430 3 mg/kg single intravenous bolus in 50% Solutol HS15 / 50%",
      "polyethylene glycol 400; plasma sampled to 10 h, liver excised at",
      "10 h."
    ),
    regions = "Republic of Korea (Kyung Hee University)",
    notes = paste(
      "n = 6 rats per group for the intravenous study (Table 1). The model",
      "itself was not fitted to individual animals: the hepatic inputs",
      "come from hepatocytes isolated from repeated-dosing rats, scaled by IVIVE",
      "Method 2 with 108 x 10^6 cells/g liver and 36 g liver/kg (Section",
      "2.5); the systemic volume and the single-adjustment rate constants",
      "are the control-group fitted values (Table 3, 'Equal to Control').",
      "Observed non-compartmental values (Table 1): AUC0-inf 1420 +/- 446",
      "ug/mL*min, CL 137.9 +/- 43.8 mL/h/kg, Vdss 378.4 +/- 119.3 mL/kg,",
      "Rb 4.53 +/- 0.65, fu,p 0.0081 +/- 0.002. Biliary and renal clearance",
      "of parent SCR430 were below 0.1% of total clearance."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Systemic and single-adjustment compartments (Table 3, 'Fitted from
    # observed data'; shared by both groups -- 'Equal to Control').
    # The systemic compartment carries the BLOOD concentration: Eq. A1
    # multiplies it by the hepatic blood flow Qh, and unbound uptake uses
    # fu,b. k12 / k21 are the paper's kin / kout into and out of the
    # single-adjustment compartment (Figure 7), first-order on amount.
    # ------------------------------------------------------------------
    lvc <- log(0.018)
    label("Systemic (central) blood compartment volume (L/kg)") # Table 3 Tissue volume 'Central' 0.018 L/kg, fitted from observed data
    lk12 <- log(2.7)
    label("Rate constant from systemic to single-adjustment compartment kin (1/h)") # Table 3 k_in 2.7 /h, fitted from observed data
    lk21 <- log(1.0)
    label("Rate constant from single-adjustment to systemic compartment kout (1/h)") # Table 3 k_out 1.0 /h, fitted from observed data

    # ------------------------------------------------------------------
    # Liver physiology (Table 3).
    # ------------------------------------------------------------------
    q_liver <- fixed(3.69)
    label("Hepatic blood flow Qh (L/h/kg)") # Table 3 Blood flow Liver 3.69 L/h/kg [ref 28]; Section 2.7 Qh 3.69 L/h/kg
    v_is_liver <- fixed(0.01)
    label("Liver extracellular space volume Vhe, sum of five units (L/kg)") # Table 3 'Extracellular space in the liver' 0.01 L/kg [ref 44]
    v_int_liver <- fixed(0.027)
    label("Hepatocyte volume Vhc, sum of five units (L/kg)") # Table 3 'Hepatocytes' 0.027 L/kg [ref 44]; with Vhe sums to Table 3 Liver 0.037 L/kg

    # ------------------------------------------------------------------
    # Blood binding (Table 1, repeated-dosing group). Section 2.9: 'The fu,b was
    # calculated from the fu,p and the Rb'.
    # ------------------------------------------------------------------
    fup <- fixed(0.0081)
    label("Fraction unbound in plasma fu,p (unitless)") # Table 1 repeated-dosing fu,p 0.0081 +/- 0.002
    bpr <- fixed(4.53)
    label("Blood-to-plasma concentration ratio Rb (unitless)") # Table 1 repeated-dosing Rb 4.53 +/- 0.65; Section 3.1

    # ------------------------------------------------------------------
    # Hepatic IVIVE inputs (Table 3 'Monte Carlo simulation parameters',
    # repeated-dosing column, reported as (mean, SD) of a lognormal distribution).
    # Section 2.11 Eqs. 17-19 parameterise that distribution by its
    # ARITHMETIC mean mu_x and SD sigma_x:
    #   mu_w = ln(mu_x^2 / sqrt(sigma_x^2 + mu_x^2)),
    #   sigma_w^2 = ln(1 + sigma_x^2 / mu_x^2).
    # The typical values below are exp(mu_w), the MEDIAN of that
    # distribution, so that exp(l<param> + eta<param>) reproduces it
    # exactly; the eta variances are sigma_w^2.
    # ------------------------------------------------------------------
    lps_inf_act <- fixed(log(9.2) - 0.5 * log(1 + (0.6 / 9.2)^2))
    label("Active hepatic uptake clearance PSinf,act, median before scaling (L/h/kg)") # Table 3 PSinf,act (mean 9.2, SD 0.6) L/h/kg; Eq. 17
    lps_dif <- fixed(log(0.86) - 0.5 * log(1 + (0.4 / 0.86)^2))
    label("Passive diffusion clearance PSdiff, median before scaling (L/h/kg)") # Table 3 PSdiff (mean 0.86, SD 0.4) L/h/kg; Eq. 17
    lcl_int_met <- fixed(log(12.3) - 0.5 * log(1 + (0.08 / 12.3)^2))
    label("Intrinsic metabolic clearance CLint,met, median (L/h/kg)") # Table 3 CLint,met (mean 12.3, SD 0.08) L/h/kg; Eq. 17
    lfu_hepa <- fixed(log(0.01) - 0.5 * log(1 + (0.008 / 0.01)^2))
    label("Fraction unbound in hepatocytes fu,hepa, median (unitless)") # Table 3 fu,hepa (mean 0.01, SD 0.008); Eq. 17; Table 1 repeated-dosing fu,hepa 0.01 +/- 0.008
    sf_ps <- fixed(1.6)
    label("Empirical scaling factor on PSinf,act and PSdiff (unitless)") # Section 3.7: 'the empirical scaling factor (1.6) was multiplied to match the in vivo total CL'

    # Lognormal Monte Carlo variability (Section 2.11, Eq. 18). The source
    # assigns it to the in vitro hepatocyte assay variability; it is the
    # only variability the paper reports.
    etalps_inf_act ~ fixed(log(1 + (0.6 / 9.2)^2))
    etalps_dif ~ fixed(log(1 + (0.4 / 0.86)^2))
    etalcl_int_met ~ fixed(log(1 + (0.08 / 12.3)^2))
    etalfu_hepa ~ fixed(log(1 + (0.008 / 0.01)^2))

    # No residual-error model is reported (the source compares Monte Carlo
    # percentiles with group means); fixed to zero rather than invented.
    propSd <- fixed(0)
    label("Proportional residual error on plasma SCR430 (fraction; not reported in source)")
    propSd_Cliver <- fixed(0)
    label("Proportional residual error on hepatocellular SCR430 (fraction; not reported in source)")
  })

  model({
    # 1. Individual parameters
    vc <- exp(lvc)
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    ps_inf_act <- exp(lps_inf_act + etalps_inf_act)
    ps_dif <- exp(lps_dif + etalps_dif)
    cl_int_met <- exp(lcl_int_met + etalcl_int_met)
    fu_hepa <- exp(lfu_hepa + etalfu_hepa)

    # 2. Derived hepatic terms. Section 3.7 applies the 1.6-fold factor to
    # both PSinf,act and PSdiff; Table 3 lists them before scaling
    # (9.2 + 0.86 = 10.06 L/h/kg reproduces the unscaled Method 2 PSinf of
    # 168.0 mL/min/kg in Table 2). Uptake (influx) is active + passive,
    # efflux is passive only (Eqs. A1-A3).
    fub <- fup / bpr
    ps_upt <- sf_ps * (ps_inf_act + ps_dif)
    ps_eff <- sf_ps * ps_dif
    v_is_seg <- v_is_liver / 5
    v_int_seg <- v_int_liver / 5

    # 3. Concentrations (mg/L = ug/mL)
    c_sys <- central / vc
    c_he1 <- is_liver1 / v_is_seg
    c_he2 <- is_liver2 / v_is_seg
    c_he3 <- is_liver3 / v_is_seg
    c_he4 <- is_liver4 / v_is_seg
    c_he5 <- is_liver5 / v_is_seg
    c_hc1 <- int_liver1 / v_int_seg
    c_hc2 <- int_liver2 / v_int_seg
    c_hc3 <- int_liver3 / v_int_seg
    c_hc4 <- int_liver4 / v_int_seg
    c_hc5 <- int_liver5 / v_int_seg

    # 4. Systemic and single-adjustment compartments (Figure 7). The
    # source prints no equation for these two; the structure is read from
    # the Figure 7 schematic: hepatic blood leaves the systemic compartment
    # into Inlet 1 and returns from Inlet 5, and kin / kout link the
    # single-adjustment compartment.
    d/dt(central) <- q_liver * (c_he5 - c_sys) - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # 5. Liver extracellular units (Eqs. A1 and A2, written on amounts:
    # d/dt(is_liver<i>) = (Vhe/5) * dChe,i/dt).
    d/dt(is_liver1) <- q_liver * (c_sys - c_he1) - ps_upt / 5 * fub * c_he1 + ps_eff / 5 * fu_hepa * c_hc1
    d/dt(is_liver2) <- q_liver * (c_he1 - c_he2) - ps_upt / 5 * fub * c_he2 + ps_eff / 5 * fu_hepa * c_hc2
    d/dt(is_liver3) <- q_liver * (c_he2 - c_he3) - ps_upt / 5 * fub * c_he3 + ps_eff / 5 * fu_hepa * c_hc3
    d/dt(is_liver4) <- q_liver * (c_he3 - c_he4) - ps_upt / 5 * fub * c_he4 + ps_eff / 5 * fu_hepa * c_hc4
    d/dt(is_liver5) <- q_liver * (c_he4 - c_he5) - ps_upt / 5 * fub * c_he5 + ps_eff / 5 * fu_hepa * c_hc5

    # 6. Hepatocellular units (Eq. A3, on amounts).
    d/dt(int_liver1) <- ps_upt / 5 * fub * c_he1 - (ps_eff + cl_int_met) / 5 * fu_hepa * c_hc1
    d/dt(int_liver2) <- ps_upt / 5 * fub * c_he2 - (ps_eff + cl_int_met) / 5 * fu_hepa * c_hc2
    d/dt(int_liver3) <- ps_upt / 5 * fub * c_he3 - (ps_eff + cl_int_met) / 5 * fu_hepa * c_hc3
    d/dt(int_liver4) <- ps_upt / 5 * fub * c_he4 - (ps_eff + cl_int_met) / 5 * fu_hepa * c_hc4
    d/dt(int_liver5) <- ps_upt / 5 * fub * c_he5 - (ps_eff + cl_int_met) / 5 * fu_hepa * c_hc5

    # 7. Observations. Plasma = systemic blood / Rb. The liver output is the
    # mean hepatocellular concentration: the measured liver concentrations
    # were corrected for the blood remaining in the tissue (Section 2.3,
    # Eq. 1), and the hepatocellular mean reproduces the Figure 8c/8d
    # median curves where the extracellular-inclusive average runs about
    # twofold above them.
    Cc <- c_sys / bpr
    Cliver <- (int_liver1 + int_liver2 + int_liver3 + int_liver4 + int_liver5) / v_int_liver

    Cc ~ prop(propSd)
    Cliver ~ prop(propSd_Cliver)
  })
}
