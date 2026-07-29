Zhou_2016_warfarin_vk2 <- function() {
  description <- "Two-drug population PK/PD model for warfarin and intravenous vitamin K2 (menatetrenone) in Japanese adults with atrial fibrillation undergoing catheter ablation. Warfarin and vitamin K2 each have a 1-compartment PK with fixed volumes-of-distribution (Vd1 = 0.183 L/kg for warfarin from Sato 2006; Vd3 = 0.051 L/kg for vitamin K2 from the Eisai product information) and fixed warfarin elimination rate (k10 = 0.0129 1/h); only the vitamin K2 elimination rate (k30) and the indirect-response PD parameters (ks, kd, IC50, Emax, EC50) were estimated from 579 INR observations in 100 patients. Warfarin inhibits clotting-factor synthesis (Emax = 1 - Cp1/(Cp1 + IC50)) while vitamin K2 stimulates it (1 + Emax_vk2 * Cp3/(Cp3 + EC50)); a binary renal-impairment indicator (CREAT >= 1.1 mg/dL in men or >= 0.8 mg/dL in women) reduces IC50 to 61.4% of normal. The model predicts thrombotest (TT, %); INR is recovered from TT via the Gogstad 1986 quadratic conversion (Equation 4)."
  reference   <- "Zhou Z, Yano I, Odaka S, Morita Y, Shizuta S, Hayano M, Kimura T, Akaike A, Inui K-i, Matsubara K. Effect of vitamin K2 on the anticoagulant activity of warfarin during the perioperative period of catheter ablation: Population analysis of retrospective clinical data. J Pharm Health Care Sci. 2016;2:17. doi:10.1186/s40780-016-0053-8. Fixed warfarin PK from Sato 2006 Jpn J Ther Drug Monit 23:10-16; vitamin K2 Vd from Eisai product information. INR <-> TT conversion from Gogstad 1986 Thromb Haemost 56:178-182."
  vignette    <- "Zhou_2016_warfarin_vk2"
  paper_specific_compartments <- c("central_vk2")

  units       <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight (baseline)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used to scale warfarin and vitamin K2 volumes-of-distribution and clearances. Median 63.8 kg, range 34.9-92.6 kg in Zhou 2016 Table 1.",
      source_name        = "WT"
    ),
    CREAT = list(
      description        = "Serum creatinine concentration (baseline)",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used together with SEXF to derive the binary renal-impairment indicator (RF) inside model(): RF = 1 when CREAT >= 1.1 mg/dL in men (SEXF = 0) or CREAT >= 0.8 mg/dL in women (SEXF = 1), otherwise 0. See Zhou 2016 Methods (RF definition paragraph) and Equation 8. Cohort median 0.8 mg/dL, range 0.5-9.6 mg/dL (Zhou 2016 Table 1). 22 of 100 patients had RF = 1.",
      source_name        = "Serum creatinine"
    ),
    SEXF = list(
      description        = "Biological sex (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Sets the sex-specific creatinine threshold for the binary renal-impairment indicator (>=1.1 mg/dL men, >=0.8 mg/dL women; Zhou 2016 Methods). Cohort distribution 70 men / 30 women (Zhou 2016 Table 1).",
      source_name        = "Sex"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 100L,
    n_observations  = 579L,
    n_studies       = 1L,
    age_range       = "31-80 years (median 64)",
    weight_range    = "34.9-92.6 kg (median 63.8)",
    sex_female_pct  = 30,
    race_ethnicity  = c(Asian = 100),
    disease_state   = "Adult Japanese patients with atrial fibrillation undergoing catheter ablation; on chronic warfarin anticoagulation withdrawn perioperatively and antagonized with intravenous vitamin K2.",
    dose_range      = "Warfarin oral maintenance 1-7 mg/day (median 3); post-operative loading 1-9 mg (median 5). Intravenous vitamin K2 total 20-70 mg (median 40) preoperatively in 76 of 100 patients (20 mg in 19, 30 mg in 2, 40 mg in 35, 60 mg in 19, 70 mg in 1).",
    regions         = "Japan (Kyoto University Hospital, January-December 2008).",
    renal_function  = "22 of 100 patients had a serum creatinine above the in-hospital reference (>=1.1 mg/dL men or >=0.8 mg/dL women); 26 patients had eGFR 30-60 mL/min/1.73 m^2 and 2 had eGFR < 30 mL/min/1.73 m^2.",
    hepatic_function = "4 of 100 patients had total bilirubin above the in-hospital reference; none were substantially elevated. 8 of 100 patients had serum albumin below the reference. Hepatic-impairment effects were not retained in the final model.",
    notes           = "Retrospective single-center cohort. Cohort and dosing details from Zhou 2016 Table 1. INR values in the analysis window were 1.0-3.0 (inclusion criterion); initial INR median 1.76, range 1.03-2.64. Concomitant amiodarone (n=4) and bucolome (n=1) were tested as CYP2C9 inhibitors on warfarin k10 but did not reach significance (-2LLD = 7.61 < 7.88)."
  )

  ini({
    # ===================================================================
    # Warfarin PK -- FIXED from Sato 2006 (the upstream Japanese popPK
    # paper). The Zhou 2016 analysis did not estimate warfarin PK; only
    # INR / TT data were available.
    #   Sato 2006: k10 = 0.0129 1/h, Vd1 = 0.183 L/kg.
    # Reparameterised as CL = k10 * Vd1 (L/(kg*h)) and Vc = Vd1 (L/kg)
    # per the lcl/lvc canonical convention; model() multiplies by WT.
    #   lcl = log(0.0129 * 0.183) = log(0.002361)
    #   lvc = log(0.183)
    # ===================================================================
    lcl <- fixed(log(0.0129 * 0.183)); label("Warfarin apparent clearance per kg body weight (L/(kg*h))")  # Sato 2006, fixed in Zhou 2016 Methods
    lvc <- fixed(log(0.183));          label("Warfarin apparent volume of distribution per kg body weight (L/kg)")  # Sato 2006, fixed in Zhou 2016 Methods

    # ===================================================================
    # Vitamin K2 PK -- Vd3 fixed from Eisai product information; k30 was
    # the only PK parameter estimated in Zhou 2016 (Table 2).
    #   Vd3 = 0.051 L/kg (fixed); k30 = 0.0194 1/h (estimated, RSE 19.2%).
    # Reparameterised analogously:
    #   lcl_vk2 = log(k30 * Vd3) = log(0.0194 * 0.051) = log(0.0009894)
    #   lvc_vk2 = fixed(log(0.051))
    # IIV on k30 (omega^2 = 0.171, CV 41.4%) maps 1:1 onto IIV on lcl_vk2
    # because lvc_vk2 is fixed.
    # ===================================================================
    lcl_vk2 <- log(0.0194 * 0.051); label("Vitamin K2 apparent clearance per kg body weight (L/(kg*h))")  # Zhou 2016 Table 2 (k30 = 0.0194 1/h; combined with Eisai Vd3 = 0.051 L/kg)
    lvc_vk2 <- fixed(log(0.051));   label("Vitamin K2 apparent volume of distribution per kg body weight (L/kg)")  # Eisai product information, fixed in Zhou 2016 Methods

    # ===================================================================
    # Indirect-response PD parameters -- Zhou 2016 Table 2 (final model).
    # State: clotting-factor activity in thrombotest-% units (TT).
    # Drug-free steady-state TT = ksyn / kd.
    # ===================================================================
    lksyn <- log(3.97);    label("Zero-order clotting-factor synthesis rate (%/h)")  # Zhou 2016 Table 2 (ks = 3.97 %/h, RSE 17.5%)
    lkd   <- log(0.0611);  label("First-order clotting-factor degradation rate (1/h)")  # Zhou 2016 Table 2 (kd = 0.0611 1/h, RSE 9.90%)
    lic50 <- log(0.604);   label("Warfarin 50% inhibitory concentration on synthesis (ug/mL)")  # Zhou 2016 Table 2 (IC50 = 0.604 ug/mL, RSE 24.5%)
    lemax <- log(0.324);   label("Vitamin K2 maximum stimulatory effect on synthesis (unitless)")  # Zhou 2016 Table 2 (Emax = 0.324, RSE 15.9%)
    lec50 <- log(5.30);    label("Vitamin K2 50% effective concentration on synthesis (ug/mL)")  # Zhou 2016 Table 2 (EC50 = 5.30 ug/mL, RSE 17.6%)

    # ===================================================================
    # Covariate effect on warfarin IC50 -- binary renal-impairment
    # indicator (Zhou 2016 Methods Equation 8 / Table 2).
    #   IC50_RF = IC50 * theta^RF, theta = 0.614 (RSE 13.9%).
    # Encoded as a log-multiplicative effect:
    #   e_creat_ic50 = log(0.614) = -0.4877
    # When the derived RF flag is 1 the factor is exp(e_creat_ic50) =
    # 0.614 (warfarin is more potent under renal impairment).
    # ===================================================================
    e_creat_ic50 <- log(0.614); label("Log-multiplicative effect of binary renal impairment on warfarin IC50 (unitless)")  # Zhou 2016 Table 2 (theta = 0.614, RSE 13.9%)

    # ===================================================================
    # Inter-individual variability -- Zhou 2016 Table 2.
    # Only ks, IC50, and k30 carry IIV (minimum-AIC eta selection).
    # omega^2 reported directly; the paper's CV% convention is
    # CV = sqrt(omega^2) * 100, not the log-normal CV formula.
    #   ks:   omega^2 = 0.0704  (CV 26.5%, RSE 25.6%)
    #   IC50: omega^2 = 0.144   (CV 37.9%, RSE 43.3%)
    #   k30:  omega^2 = 0.171   (CV 41.4%, RSE 85.3%)
    # ===================================================================
    etalksyn   ~ 0.0704
    etalic50   ~ 0.144
    etalcl_vk2 ~ 0.171

    # ===================================================================
    # Residual error -- Zhou 2016 Methods Equation 6:
    #   TT_obs = TT_pred * exp(eps), eps ~ N(0, sigma^2).
    # This is a log-normal residual error on TT. Table 2:
    #   sigma^2 = 0.0798 (CV 28.2%, RSE 11.8%); expSd = sqrt(0.0798).
    # ===================================================================
    expSd <- sqrt(0.0798); label("Log-normal residual error SD on thrombotest (unitless, on log scale)")  # Zhou 2016 Table 2 (sigma^2 = 0.0798, CV 28.2%, RSE 11.8%)
  })

  model({
    # -----------------------------------------------------------------
    # Derived covariate: binary renal-impairment indicator RF, per
    # Zhou 2016 Methods. RF = 1 when serum creatinine exceeds the
    # sex-specific in-hospital reference value (>=1.1 mg/dL in men,
    # >=0.8 mg/dL in women), otherwise 0.
    # -----------------------------------------------------------------
    rf <- SEXF * (CREAT >= 0.8) + (1 - SEXF) * (CREAT >= 1.1)

    # -----------------------------------------------------------------
    # Individual PK parameters. Both warfarin (Vd1) and vitamin K2
    # (Vd3) are per-kg volumes from the source paper, so multiply by
    # WT to obtain L. Warfarin clearance and volume are fixed; only
    # vitamin K2 clearance has IIV (eta on k30 -> eta on cl_vk2 with
    # vc_vk2 fixed).
    # -----------------------------------------------------------------
    cl     <- exp(lcl)               * WT
    vc     <- exp(lvc)               * WT
    cl_vk2 <- exp(lcl_vk2 + etalcl_vk2) * WT
    vc_vk2 <- exp(lvc_vk2)           * WT

    # -----------------------------------------------------------------
    # Individual PD parameters.
    # -----------------------------------------------------------------
    ksyn  <- exp(lksyn + etalksyn)
    kd    <- exp(lkd)
    ic50  <- exp(lic50 + etalic50) * exp(e_creat_ic50 * rf)
    emax  <- exp(lemax)
    ec50  <- exp(lec50)

    # -----------------------------------------------------------------
    # ODE system. Both drugs use 1-compartment PK with instantaneous
    # bolus input into the central compartment (warfarin orally,
    # vitamin K2 intravenously -- the paper treats absorption as
    # instantaneous; Methods Equations 1-2). Plasma concentrations
    # Cp1 (warfarin) and Cp3 (vitamin K2) drive the indirect-response
    # clotting-factor model.
    #
    #   d/dt(central)     = -kel    * central       (warfarin)
    #   d/dt(central_vk2) = -kel_vk2 * central_vk2  (vitamin K2)
    #   d/dt(effect)      = ksyn * (1 - Cp1/(Cp1+IC50) + Emax*Cp3/(Cp3+EC50))
    #                     - kd * effect
    #
    # Zhou 2016 Methods Equation 9 (covariate-adjusted indirect-response
    # equation with IC50 -> IC50 * theta^RF).
    # -----------------------------------------------------------------
    d/dt(central)     <- -(cl     / vc)     * central
    d/dt(central_vk2) <- -(cl_vk2 / vc_vk2) * central_vk2

    Cp1 <- central     / vc
    Cp3 <- central_vk2 / vc_vk2

    inhibition  <- 1 - Cp1 / (Cp1 + ic50)
    stimulation <- emax * Cp3 / (Cp3 + ec50)

    d/dt(effect) <- ksyn * (inhibition + stimulation) - kd * effect

    # Drug-free steady-state baseline TT = ksyn / kd. The simulation
    # warms up on chronic warfarin maintenance dosing to reach the
    # treated steady state before any perioperative perturbation.
    effect(0) <- ksyn / kd

    # -----------------------------------------------------------------
    # Observation: thrombotest (TT, %). INR is recovered from TT in
    # the vignette by inverting Zhou 2016 Methods Equation 4
    # (Gogstad 1986 calibration), since INR is the clinical readout
    # but the model is parameterised on TT.
    # -----------------------------------------------------------------
    TT <- effect
    TT ~ lnorm(expSd)
  })
}
