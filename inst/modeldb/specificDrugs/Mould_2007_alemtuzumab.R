Mould_2007_alemtuzumab <- function() {
  description <- "Two-compartment population PK model with Michaelis-Menten elimination for alemtuzumab in B-cell chronic lymphocytic leukaemia (Mould 2007)"
  reference <- "Mould DR, Baumann A, Kuhlmann J, Keating MJ, Weitman S, Hillmen P, Brettman LR, Reif S, Bonate PL. Population pharmacokinetics-pharmacodynamics of alemtuzumab (Campath) in patients with chronic lymphocytic leukaemia and its link to treatment response. Br J Clin Pharmacol. 2007;64(3):278-291. doi:10.1111/j.1365-2125.2007.02914.x"
  vignette <- "Mould_2007_alemtuzumab"
  units <- list(time = "hour", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WBC = list(
      description        = "Total white blood cell count (peripheral, time-varying biomarker of B-CLL tumour burden)",
      units              = "10^9 cells/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Power covariate on Vmax with reference 10 x 10^9/L. In the Mould 2007 cohort baseline WBC was markedly elevated (median 37.8 x 10^9/L, range 1.3-522) because circulating leukaemic B-cells dominate the count, and WBC falls during alemtuzumab treatment as the clone is depleted; simulations must therefore supply WBC at every observation time.",
      source_name        = "WBC"
    )
  )

  population <- list(
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
    baseline_wbc     = "median 37.8 x 10^9/L (range 1.3-522) — reflects B-CLL leukaemic burden",
    notes            = "Baseline demographics per Mould 2007 Table 1 (four pooled studies: CAM002, CAM005, CAM211, CAM213). Weight was tested but not retained as a covariate; WBC was the only covariate retained in the final model."
  )

  ini({
    # Structural parameters (Mould 2007 Table 2, final estimates).
    # Paper reports Vmax in ug/h and Km in ug/L; values kept in those units
    # inside ini() and converted to (mg/h) and (mg/L) inside model() so that
    # mg dosing gives concentrations in mg/L (= ug/mL).
    lvmax <- log(1020); label("Typical Vmax at reference WBC = 10 x 10^9/L (ug/h)")  # Mould 2007 Table 2
    lkm   <- log(338);  label("Michaelis-Menten constant Km (ug/L)")                  # Mould 2007 Table 2
    lvc   <- log(11.3); label("Central volume of distribution V1 (L)")                # Mould 2007 Table 2
    lq    <- log(1.05); label("Intercompartmental clearance Q (L/h)")                 # Mould 2007 Table 2
    lvp   <- log(41.5); label("Peripheral volume of distribution V2 (L)")             # Mould 2007 Table 2

    # Covariate effect: Vmax = TVVmax * (WBC / 10)^e_wbc_vmax (power form)
    e_wbc_vmax <- 0.194; label("Power exponent for WBC on Vmax (unitless)")           # Mould 2007 Table 2

    # Inter-individual variability (simple diagonal Omega, log-normal).
    # Paper reports ISV as %CV; omega^2 = log(CV^2 + 1).
    etalvmax ~ 0.097460  # 32%  CV on Vmax (Mould 2007 Table 2)
    etalkm   ~ 1.132194  # 145% CV on Km   (Mould 2007 Table 2)
    etalvc   ~ 0.534032  # 84%  CV on V1   (Mould 2007 Table 2)
    etalvp   ~ 1.436193  # 179% CV on V2   (Mould 2007 Table 2)

    # Combined residual error (Mould 2007 Table 2).
    # CCV = 37.2% proportional; additive 64.73 ug/L = 0.06473 ug/mL.
    propSd <- 0.372;   label("Proportional residual error (fraction)")                # Mould 2007 Table 2
    addSd  <- 0.06473; label("Additive residual error (ug/mL)")                       # Mould 2007 Table 2 (64.73 ug/L / 1000)
  })
  model({
    # Individual PK parameters. Vmax has a power covariate on WBC
    # (reference 10 x 10^9/L) and is converted from ug/h to mg/h so that
    # mg-denominated dosing yields mg/L (= ug/mL) concentrations.
    vmax <- exp(lvmax + etalvmax) * (WBC / 10)^e_wbc_vmax / 1000
    km   <- exp(lkm   + etalkm)   / 1000
    vc   <- exp(lvc   + etalvc)
    vp   <- exp(lvp   + etalvp)
    q    <- exp(lq)

    k12 <- q / vc
    k21 <- q / vp

    Cc_central <- central / vc

    # Two-compartment IV input model with Michaelis-Menten elimination
    # from the central compartment.
    d/dt(central)     <- -vmax * Cc_central / (km + Cc_central) - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    Cc <- Cc_central
    Cc ~ add(addSd) + prop(propSd)
  })
}
