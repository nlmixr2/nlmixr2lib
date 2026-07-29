Ye_2017_ethaselen <- function() {
  description <- paste(
    "Preclinical (mouse, BALB/c nude with A549 NSCLC xenograft).",
    "Integrated dose-biomarker-response PD model for the thioredoxin",
    "reductase (TrxR) inhibitor ethaselen (Ye et al. 2017). The TrxR",
    "biomarker is described by an indirect-response (IDR) turnover",
    "in which the zero-order production Kin is linearly amplified by",
    "the instantaneous natural tumor growth rate (linear correction",
    "factor gamma1) and the first-order degradation Kout is increased",
    "by a sigmoidal Emax function of the current administered ethaselen",
    "dose (Smax, SC50, Hill = gamma2). Tumor volume follows a smooth",
    "exponential-to-linear growth law (paper Eq 5:",
    "dX/dt = 2*lambda0*lambda1*X / (lambda1 + 2*lambda0*X)) tempered",
    "by a zero-order Emax killing rate driven by the TrxR-inhibition",
    "ratio P = 1 - TrxR_treatment / TrxR_control (paper Eq 7). The",
    "control TrxR trajectory is carried internally as a shadow state",
    "(trxr_ctrl) so P is defined per-subject without requiring an",
    "external control-arm simulation. No pharmacokinetic compartment",
    "is included; the paper acknowledges ethaselen plasma",
    "concentrations were not measured. The current daily dose enters",
    "the model through the time-varying covariate DOSE (mg/kg/day),",
    "which the published study toggles between 0 (vehicle /",
    "off-treatment) and one of {36, 72, 108} mg/kg/day for days 0-9",
    "(oral gavage QD x10 d)."
  )
  reference <- paste(
    "Ye SF, Li J, Ji SM, Zeng HH, Lu W.",
    "Dose-biomarker-response modeling of the anticancer effect of",
    "ethaselen in a human non-small cell lung cancer xenograft mouse",
    "model. Acta Pharmacologica Sinica 2017;38(2):223-232.",
    "doi:10.1038/aps.2016.114.",
    "Published online 2016-12-05; print issue 2017.",
    sep = " "
  )
  vignette <- "Ye_2017_ethaselen"
  paper_specific_compartments <- c("trxr", "trxr_ctrl", "tumor_volume")

  units <- list(
    time          = "day",
    dosing        = "mg/kg (administered as a per-record DOSE covariate; no PK ODE)",
    concentration = "n/a (TrxR activity in U/mL and tumor volume in mm^3 are the two observed outputs)"
  )

  covariateData <- list(
    DOSE = list(
      description        = "Time-varying current administered daily dose of ethaselen (mg/kg/day). Drives the sigmoidal Emax inhibition of TrxR Kout. Set to 0 outside the dosing window (vehicle control and post-treatment days).",
      units              = "mg/kg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Published study uses DOSE in {0, 36, 72, 108} mg/kg/day during days 0-9 (oral gavage QD x10 d) and DOSE = 0 thereafter. The drug-effect term collapses to zero when DOSE = 0 because DOSE^hill = 0 (hill = 2.29 > 0), so the shadow control state trxr_ctrl coincides with the treated state trxr in vehicle subjects and P = 0.",
      source_name        = "D"
    )
  )

  population <- list(
    species        = "mouse (female BALB/c nude, 5 weeks old, with subcutaneous A549 human non-small cell lung cancer xenograft)",
    n_subjects     = 160L,
    n_studies      = 2L,
    age_range      = "5 weeks at randomisation",
    weight_range   = "not reported",
    sex_female_pct = 100,
    race_ethnicity = NA,
    disease_state  = "subcutaneous A549 NSCLC xenograft (5x10^6 cells / 0.2 mL inoculated into the right armpit); modelling started when tumor volume reached approximately 100 mm^3",
    dose_range     = "vehicle (0.5% CMC-Na, pH 7.4), 36, 72, or 108 mg/kg ethaselen ig QD x10 d",
    regions        = "preclinical (in-vivo xenograft); Peking University Health Science Center animal facility",
    notes          = "Two sub-studies pooled into one model dataset: (i) TrxR biomarker assay -- 132 mice (4 groups x 33 mice; n=3 sacrificed per day per group across 10 dosing days plus follow-up); (ii) tumor-volume study -- 28 mice (4 groups x 7 mice; volumes recorded daily by Vernier calipers, V = (length x width^2) / 2). Total = 160 mice. Tumor TrxR was measured by the DTNB-reducing assay (linear range 0-1000 U/mL). The integrated dose-biomarker-response model in Table 1 was fit to all four dose arms using NONMEM 7.1.2 / FOCEI with PsN 3.5.3 and Pirana 2.8.0; validation via 1000-replicate VPC."
  )

  ini({
    # ---- TrxR turnover (indirect-response) -------------------------------
    # Table 1: Kin = 8.27 U/mL/day, RSE 42.4%, IIV (CV) = 12.4%.
    lkin    <- log(8.27)
    label("TrxR zero-order production rate Kin (U/mL/day)")              # Ye 2017 Table 1, Kin = 8.27 (RSE 42.4%)
    # Table 1: baseline TrxR activity Base = 39.7 U/mL, RSE 8.6%, IIV (CV) = 9.6%.
    # Kout is derived inside model() as Kin/rbase to maintain steady state.
    lrbase  <- log(39.7)
    label("Baseline TrxR activity (U/mL)")                                # Ye 2017 Table 1, Base = 39.7 (RSE 8.60%)

    # ---- Linear correction of Kin by tumor growth rate -------------------
    # Table 1 reports gamma1 = 0.021 with units "d/mm" and RSE 16.5%.
    # Applied multiplicatively to Kin per the paper's prose
    # ("the linear correction factor gamma1"). The published unit
    # "d/mm" appears to be a typographical compression of "d/mm^3";
    # the value is reproduced exactly as published. See vignette
    # Errata for the unit reconciliation.
    lgamma  <- log(0.021)
    label("Linear correction factor on Kin from natural tumor growth rate gamma1 (d/mm^3)")  # Ye 2017 Table 1, gamma1 = 0.021 (RSE 16.5%)

    # ---- Sigmoidal Emax effect of dose on Kout (TrxR degradation) --------
    # Table 1: Smax = 5.95 (RSE 31.9%), SC50 = 136 mg/kg (RSE 25.2%),
    # gamma2 = Hill coefficient = 2.29 (RSE 17.3%). No IIV on any of
    # the three (NE in Table 1). The drug effect multiplies Kout:
    # Kout_eff = Kout * (1 + Smax * DOSE^gamma2 / (SC50^gamma2 +
    # DOSE^gamma2)).
    lsmax   <- log(5.95)
    label("Maximal multiplicative TrxR-Kout stimulation Smax (unitless)")  # Ye 2017 Table 1, Smax = 5.95 (RSE 31.9%)
    lsc50   <- log(136)
    label("Dose of ethaselen at half-maximal TrxR-Kout stimulation SC50 (mg/kg)")  # Ye 2017 Table 1, SC50 = 136 (RSE 25.2%)
    lhill   <- log(2.29)
    label("Hill coefficient on dose-driven TrxR-Kout stimulation gamma2 (unitless)")  # Ye 2017 Table 1, gamma2 = 2.29 (RSE 17.3%)

    # ---- Exponential-to-linear tumor natural growth (paper Eq 5) ---------
    # dX/dt = 2 * lambda0 * lambda1 * X / (lambda1 + 2 * lambda0 * X).
    # Behaviour: at small X the rate approaches 2*lambda0 * X
    # (exponential, doubling time ln(2) / (2 * lambda0)); at large X
    # the rate saturates at lambda1 (linear). The transition occurs at
    # X = lambda1 / (2 * lambda0).
    # Table 1: lambda0 = 0.704/day (RSE 31.4%, no IIV);
    # lambda1 = 321 mm^3/day (RSE 11.5%, IIV CV 32.7%).
    llambda0 <- log(0.704)
    label("Exponential growth rate constant lambda0 (1/day)")             # Ye 2017 Table 1, lambda0 = 0.704 (RSE 31.4%)
    llambda1 <- log(321)
    label("Linear growth rate constant lambda1 (mm^3/day)")               # Ye 2017 Table 1, lambda1 = 321 (RSE 11.5%)

    # ---- Initial tumor volume --------------------------------------------
    # Table 1: W = 103 mm^3, RSE 3.9%, IIV (CV) = 14.2%.
    lrbase_tumor <- log(103)
    label("Initial tumor volume W (mm^3)")                                # Ye 2017 Table 1, W = 103 (RSE 3.9%)

    # ---- Drug-effect tumor eradication (paper Eq 8) ----------------------
    # Killing rate = Emax * P / (EC50 + P) where P is the relative TrxR
    # inhibition ratio P = 1 - TrxR_treatment / TrxR_control (paper
    # Eq 7). Simple Emax with no Hill on P (none reported in Table 1).
    # Table 1: Emax = 130 mm^3/day (RSE 4.8%, no IIV);
    # EC50 = 0.0676 ratio (RSE 23.1%, IIV CV 19.3%).
    lemax  <- log(130)
    label("Maximal tumor eradication rate Emax (mm^3/day)")               # Ye 2017 Table 1, Emax = 130 (RSE 4.8%)
    lec50  <- log(0.0676)
    label("TrxR-inhibition ratio P at half-maximal tumor eradication EC50 (ratio)")  # Ye 2017 Table 1, EC50 = 0.0676 (RSE 23.1%)

    # ---- Inter-individual variability ------------------------------------
    # Table 1 reports IIV as CV%. Variances on the log scale are
    # omega^2 = log(CV^2 + 1).
    etalkin         ~ 0.01526  # log(0.124^2 + 1); CV = 12.4% (Kin)
    etalrbase       ~ 0.00917  # log(0.096^2 + 1); CV = 9.6%  (Base)
    etallambda1     ~ 0.10160  # log(0.327^2 + 1); CV = 32.7% (lambda1)
    etalrbase_tumor ~ 0.01996  # log(0.142^2 + 1); CV = 14.2% (W)
    etalec50        ~ 0.03657  # log(0.193^2 + 1); CV = 19.3% (EC50)

    # ---- Residual error --------------------------------------------------
    # Table 1 reports a single Err_pro (20.22%, RSE 13.6%) and a single
    # Err_add (141 mm^3, RSE 78%). The mm^3 units on the additive error
    # identify it as applying to tumor volume. The paper does not
    # separately report a TrxR residual-error structure for the
    # integrated model; the proportional error of 20.22% is reused
    # here for the TrxR endpoint as the only published magnitude. See
    # vignette Errata for the assumption.
    propSd_trxr         <- 0.2022
    label("Proportional residual SD on TrxR activity (fraction)")          # Ye 2017 Table 1, Err_pro = 20.22% (RSE 13.6%); applied to TrxR (assumption)
    propSd_tumor_volume <- 0.2022
    label("Proportional residual SD on tumor volume (fraction)")           # Ye 2017 Table 1, Err_pro = 20.22% (RSE 13.6%)
    addSd_tumor_volume  <- 141
    label("Additive residual SD on tumor volume (mm^3)")                   # Ye 2017 Table 1, Err_add = 141 (RSE 78%)
  })

  model({
    # ---- Individual structural parameters (lognormal IIV where present) --
    kin          <- exp(lkin   + etalkin)
    rbase        <- exp(lrbase + etalrbase)
    gamma        <- exp(lgamma)
    smax         <- exp(lsmax)
    sc50         <- exp(lsc50)
    hill         <- exp(lhill)
    lambda0      <- exp(llambda0)
    lambda1      <- exp(llambda1 + etallambda1)
    rbase_tumor  <- exp(lrbase_tumor + etalrbase_tumor)
    emax         <- exp(lemax)
    ec50         <- exp(lec50 + etalec50)

    # ---- Derived rate constants ------------------------------------------
    # Kout is constrained by steady-state condition Kin = Kout * Base, so
    # Kout = Kin / Base. The system rests at the reported baseline TrxR
    # activity when DOSE = 0 and tumor growth rate is zero.
    kout <- kin / rbase

    # ---- Natural tumor growth rate (paper Eq 5) --------------------------
    # dX/dt_natural = 2 * lambda0 * lambda1 * X / (lambda1 + 2 * lambda0 * X)
    # is the natural (undisturbed) tumor growth rate. Used both inside
    # the Kin amplification term (linear correction by tumor growth)
    # and as the production half of the treated-tumor ODE.
    growth_rate_natural <- (2.0 * lambda0 * lambda1 * tumor_volume) /
                           (lambda1 + 2.0 * lambda0 * tumor_volume)

    # ---- Kin amplified by natural tumor growth rate (paper text) ---------
    # The paper states "the total TrxR generation rates were affected by
    # the natural growth rates of the tumor" with "linear correction
    # factor gamma1". The form is multiplicative: Kin_eff = Kin * (1 +
    # gamma * growth_rate_natural). See vignette Errata for the choice
    # of multiplicative vs. additive linear correction.
    kin_eff <- kin * (1.0 + gamma * growth_rate_natural)

    # ---- Sigmoidal Emax effect of current dose on Kout -------------------
    # When DOSE = 0 the numerator DOSE^hill = 0 (hill = 2.29 > 0) so
    # drug_effect = 0; vehicle / off-treatment records collapse to the
    # untreated Kout. The 1e-30 in the denominator is defensive against
    # 0/0 at DOSE = 0 (the numerator is also 0).
    drug_effect <- smax * DOSE^hill / (sc50^hill + DOSE^hill + 1.0e-30)
    kout_eff    <- kout * (1.0 + drug_effect)

    # ---- TrxR-inhibition ratio P (paper Eq 7) ----------------------------
    # P = 1 - TrxR_treatment / TrxR_control. The control TrxR is carried
    # as the shadow state trxr_ctrl, integrated with the same kin_eff
    # but with the unmodified kout (no drug effect). This makes P
    # per-subject and well-defined under forward simulation without
    # requiring an external control-arm trajectory. Floored at 0 so
    # transient overshoots of trxr above trxr_ctrl (numerically
    # possible at simulation startup or under stochastic IIV) do not
    # drive a negative kill rate.
    p_inhib_raw <- 1.0 - trxr / trxr_ctrl
    p_inhib     <- p_inhib_raw * (p_inhib_raw > 0.0)

    # ---- Tumor eradication driven by P (paper Eq 8) ----------------------
    kill_rate <- emax * p_inhib / (ec50 + p_inhib + 1.0e-30)

    # ---- ODE system ------------------------------------------------------
    # State variables:
    #   trxr          = TrxR activity in the treated subject (U/mL)
    #   trxr_ctrl     = TrxR activity in the no-drug counterfactual
    #                   that drives the P = 1 - trxr/trxr_ctrl ratio
    #                   (U/mL; shadow state, not an observation)
    #   tumor_volume  = subcutaneous tumor volume (mm^3)
    d/dt(trxr)         <- kin_eff - kout_eff * trxr
    d/dt(trxr_ctrl)    <- kin_eff - kout     * trxr_ctrl
    d/dt(tumor_volume) <- growth_rate_natural - kill_rate

    # ---- Initial conditions ----------------------------------------------
    trxr(0)         <- rbase
    trxr_ctrl(0)    <- rbase
    tumor_volume(0) <- rbase_tumor

    # ---- Observation and residual error ----------------------------------
    trxr         ~ prop(propSd_trxr)
    tumor_volume ~ add(addSd_tumor_volume) + prop(propSd_tumor_volume)
  })
}
