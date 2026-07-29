Mould_2007_alemtuzumab_wbc <- function() {
  description <- "Coupled population PK-PD model for alemtuzumab in B-cell chronic lymphocytic leukaemia (Mould 2007): the two-compartment Michaelis-Menten PK from Mould 2007 Table 2 driven by the simulated WBC state via Vmax = TVVmax * (WBC/10)^0.194, joined to an indirect-response model on WBC (stimulation of Kout by alemtuzumab; Mould 2007 Table 3). WBC is a state variable initialised per subject at Kin/Kout."
  reference <- "Mould DR, Baumann A, Kuhlmann J, Keating MJ, Weitman S, Hillmen P, Brettman LR, Reif S, Bonate PL. Population pharmacokinetics-pharmacodynamics of alemtuzumab (Campath) in patients with chronic lymphocytic leukaemia and its link to treatment response. Br J Clin Pharmacol. 2007;64(3):278-291. doi:10.1111/j.1365-2125.2007.02914.x"
  vignette <- "Mould_2007_alemtuzumab"
  units <- list(time = "hour", dosing = "mg", concentration = "ug/mL", wbc = "10^9 cells/L")

  covariateData <- list()

  population <- list(
    species          = "human",
    n_subjects       = 67,
    n_studies        = 4,
    age_range        = "41-75 years",
    age_median       = "59 years",
    weight_range     = "45-167 kg",
    weight_median    = "72 kg",
    sex_female_pct   = 26.9,
    race_ethnicity   = "Not reported in the published analysis",
    disease_state    = "B-cell chronic lymphocytic leukaemia (B-CLL), mostly relapsed/refractory",
    dose_range       = "3-240 mg alemtuzumab as 2-h IV infusion; escalation 3 -> 10 -> 30 mg then 30 mg three times weekly (CAM005, CAM213) or 7.5-240 mg weekly for 4 weeks (CAM002)",
    regions          = "United States and Europe (pooled Schering AG clinical studies)",
    baseline_wbc     = "median 37.8 x 10^9/L (range 1.3-522); the typical-value baseline implied by the PD parameters is Kin/Kout = 1.56 / 0.029 = 53.8 x 10^9/L (Mould 2007 Results, PD section)",
    notes            = "Same pooled four-study cohort (CAM002, CAM005, CAM211, CAM213) used for the PK model in Mould_2007_alemtuzumab (PK only). The PD layer was fit sequentially: 67 patients, 1067 WBC observations (Methods, Database assembly). No covariates were tested on PD parameters (Methods, Pharmacodynamic modelling)."
  )

  ini({
    # =========================================================================
    # PK layer -- two-compartment with Michaelis-Menten elimination from
    # central. Final estimates from Mould 2007 Table 2. Paper reports Vmax in
    # ug/h and Km in ug/L; values are kept in those units inside ini() and
    # converted to (mg/h) and (mg/L) inside model() so that mg-denominated
    # dosing yields concentrations in mg/L (= ug/mL).
    # =========================================================================
    lvmax <- log(1020); label("Typical Vmax at reference WBC = 10 x 10^9/L (ug/h)")  # Mould 2007 Table 2 (TVVmax = 1020)
    lkm   <- log(338);  label("Michaelis-Menten constant Km (ug/L)")                  # Mould 2007 Table 2 (Km = 338)
    lvc   <- log(11.3); label("Central volume of distribution V1 (L)")                # Mould 2007 Table 2 (V1 = 11.3)
    lq    <- log(1.05); label("Intercompartmental clearance Q (L/h)")                 # Mould 2007 Table 2 (Q = 1.05)
    lvp   <- log(41.5); label("Peripheral volume of distribution V2 (L)")             # Mould 2007 Table 2 (V2 = 41.5)

    # Covariate effect: Vmax = TVVmax * (WBC / 10)^e_wbc_vmax (power form).
    # In this coupled model WBC is the simulated PD state (not a data column),
    # so the same coefficient and reference value (10 x 10^9/L) are reused.
    e_wbc_vmax <- 0.194; label("Power exponent for WBC on Vmax (unitless)")           # Mould 2007 Table 2 (WBC_factor = 0.194)

    # IIV on PK parameters (Mould 2007 Table 2). Paper reports ISV as %CV;
    # omega^2 = log(CV^2 + 1) for the log-normal random effects.
    etalvmax ~ 0.097460  # 32%  CV on Vmax (Mould 2007 Table 2)
    etalkm   ~ 1.132194  # 145% CV on Km   (Mould 2007 Table 2)
    etalvc   ~ 0.534032  # 84%  CV on V1   (Mould 2007 Table 2)
    etalvp   ~ 1.436193  # 179% CV on V2   (Mould 2007 Table 2)

    # PK combined residual error (Mould 2007 Table 2).
    # CCV = 37.2% proportional; additive 64.73 ug/L = 0.06473 ug/mL.
    propSd <- 0.372;   label("PK proportional residual error (fraction)")             # Mould 2007 Table 2 (CCV = 37.2%)
    addSd  <- 0.06473; label("PK additive residual error (ug/mL)")                    # Mould 2007 Table 2 (64.73 ug/L / 1000)

    # =========================================================================
    # PD layer -- indirect-response model on WBC with stimulation of Kout
    # (Mould 2007 Table 3). Final estimates from Table 3. Paper EC50 is in
    # ug/L; converted to mg/L inside model() to match the PK Cc units.
    # =========================================================================
    lemax <- log(18.2); label("Emax: maximum stimulation of WBC loss (unitless)")     # Mould 2007 Table 3 (Emax = 18.2)
    lec50 <- log(306);  label("EC50: alemtuzumab concentration at half-max effect (ug/L)") # Mould 2007 Table 3 (EC50 = 306)
    lkin  <- log(1.56); label("Kin: zero-order WBC production rate (10^9 cells/L/h)") # Mould 2007 Table 3 (Kin = 1.56)
    lkout <- log(0.029); label("Kout: first-order WBC elimination rate constant (1/h)") # Mould 2007 Table 3 (Kout = 0.029)

    # IIV on PD parameters (Mould 2007 Table 3). %CV -> omega^2 = log(CV^2 + 1).
    # Paper notes the PD model contained a correlation term between Emax and
    # EC50 but does not report the correlation value; the diagonal form is
    # used here and the missing covariance is recorded in vignette Errata.
    # Kout IIV is reported as Not estimated in Table 3 -> wrapped in fixed(0).
    etalemax ~ 1.939168    # 244% CV on Emax (Mould 2007 Table 3)
    etalec50 ~ 4.111842    # 775% CV on EC50 (Mould 2007 Table 3)
    etalkin  ~ 1.376048    # 172% CV on Kin  (Mould 2007 Table 3)
    etalkout ~ fixed(0)    # Mould 2007 Table 3 Kout IIV reported as Not estimated

    # PD residual error (Mould 2007 Table 3): additive, 15.6 x 10^9 cells/L.
    addSd_WBC <- 15.6; label("PD additive residual SD on WBC (10^9 cells/L)")         # Mould 2007 Table 3 (additive 15.6)
  })

  model({
    # ---- PD individual parameters ----
    emax <- exp(lemax + etalemax)
    # EC50 is in ug/L per Mould 2007; PK Cc below is in mg/L, so convert.
    ec50 <- exp(lec50 + etalec50) / 1000
    kin  <- exp(lkin  + etalkin)
    kout <- exp(lkout + etalkout)
    wbc_baseline <- kin / kout

    # ---- PK individual parameters ----
    # Vmax depends on the simulated WBC state (not a data covariate). Paper
    # reports Vmax in ug/h and Km in ug/L; convert to mg/h and mg/L so that
    # mg-denominated dosing yields concentrations in mg/L (= ug/mL).
    vmax <- exp(lvmax + etalvmax) * (WBC / 10)^e_wbc_vmax / 1000
    km   <- exp(lkm   + etalkm)   / 1000
    vc   <- exp(lvc   + etalvc)
    vp   <- exp(lvp   + etalvp)
    q    <- exp(lq)

    k12 <- q / vc
    k21 <- q / vp

    Cc_central <- central / vc

    # ---- ODE system ----
    # Two-compartment IV input model with Michaelis-Menten elimination from
    # the central compartment; coupled stimulatory-loss indirect response
    # model on WBC (Mould 2007 Results, PD equation): dWBC/dt = Kin - Kout *
    # (1 + Emax * C / (EC50 + C)) * WBC.
    d/dt(central)     <- -vmax * Cc_central / (km + Cc_central) - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(WBC)         <-  kin - kout * (1 + emax * Cc_central / (ec50 + Cc_central)) * WBC

    # ---- Initial condition ----
    # Steady-state baseline WBC = Kin / Kout in the absence of drug.
    WBC(0) <- wbc_baseline

    # ---- Observations and error ----
    Cc <- Cc_central
    Cc  ~ add(addSd)     + prop(propSd)
    WBC ~ add(addSd_WBC)
  })
}
