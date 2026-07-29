Schipani_2013_atazanavir_ritonavir <- function() {
  description <- paste(
    "Simultaneous one-compartment first-order-absorption popPK model for oral",
    "atazanavir (ATV) and ritonavir (RTV) in 30 HIV-infected adults receiving",
    "ATV/RTV 300/100 mg once daily, with a direct sigmoidal-Emax inhibition of",
    "ATV apparent clearance by RTV plasma concentration (Imax = 0.988,",
    "IC50 = 0.221 mg/L). Both drugs share a one-compartment structure with",
    "first-order absorption and an absorption lag time; ka values are fixed",
    "to the separate-model final estimates (ATV ka = 1.81 1/h, RTV ka =",
    "0.898 1/h) because joint estimation produced numerical instability.",
    "Inter-individual variability is carried on V/F for both drugs and on",
    "CL/F for RTV (correlated with V/F RTV, rho = 0.75); ATV CL/F is fitted",
    "without IIV. Demographic covariates and tenofovir co-administration",
    "were tested and none retained (Schipani 2013)."
  )
  reference <- paste(
    "Schipani A, Dickinson L, Boffito M, Austin R, Owen A, Back D, Khoo S,",
    "Davies G. Simultaneous Population Pharmacokinetic Modelling of Atazanavir",
    "and Ritonavir in HIV-Infected Adults and Assessment of Different Dose",
    "Reduction Strategies. J Acquir Immune Defic Syndr. 2013;62(1):60-66.",
    "doi:10.1097/QAI.0b013e3182737231."
  )
  vignette <- "Schipani_2013_atazanavir_ritonavir"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 30L,
    n_studies      = 3L,
    n_observations = 600L,
    age_range      = "22-62 years",
    age_median     = "43 years",
    weight_range   = "46-110 kg",
    weight_median  = "75.5 kg",
    sex_female_pct = 10,
    race_ethnicity = "Predominantly white; 5 Black-Africans and 3 Hispanics among 30 subjects (Schipani 2013 Methods, Patients section)",
    disease_state  = "HIV-1 infection on a regimen containing ATV/RTV 300/100 mg once daily; median baseline HIV viral load 61 copies/mL (range <50-72)",
    dose_range     = "ATV 300 mg + RTV 100 mg orally once daily (study regimen). Dose-reduction simulations in the paper additionally evaluated ATV/RTV 300/50, 200/50, and 200/100 mg once daily.",
    regions        = "United Kingdom (St Stephen's Centre, Chelsea and Westminster Foundation Trust, London)",
    notes          = paste(
      "Pooled cohort from three previously reported studies (Schipani 2013",
      "Methods refs 15-17). 27 of 30 subjects were male. 5 of 30 were",
      "receiving tenofovir 300 mg once daily; TDF co-administration was",
      "tested as a covariate on ATV CL/F and not retained. Median 11 samples",
      "per patient on a single occasion. ATV and RTV measured by HPLC-MS/MS",
      "(Schipani 2013 Methods ref 18)."
    )
  )

  ini({
    # ============================================================
    # Atazanavir (parent / substrate) structural parameters
    # ----- Schipani 2013 Table 2 -----
    # ============================================================
    lcl    <- log(16.6)
    label("Atazanavir apparent clearance in the absence of ritonavir, CL0/F (L/h)")  # Table 2: CL/F ATV = 16.6 L/h, RSE 7%
    lvc    <- log(106)
    label("Atazanavir apparent volume of distribution, V/F (L)")                     # Table 2: V/F ATV = 106 L, RSE 7%
    lka    <- fixed(log(1.81))
    label("Atazanavir first-order absorption rate, ka (1/h, FIXED)")                 # Table 2: ka ATV (per hour) fixed = 1.81 (final estimate of separate ATV model in Table 1)
    ltlag  <- log(0.87)
    label("Atazanavir absorption lag time, Tlag (h)")                                # Table 2: Lag T ATV = 0.87 h, RSE 2%

    # ============================================================
    # Ritonavir (booster / perpetrator) structural parameters
    # ----- Schipani 2013 Table 2 -----
    # ============================================================
    lcl_rtv    <- log(13.2)
    label("Ritonavir apparent clearance, CL/F (L/h)")                                # Table 2: CL/F RTV = 13.2 L/h, RSE 12%
    lvc_rtv    <- log(124)
    label("Ritonavir apparent volume of distribution, V/F (L)")                      # Table 2: V/F RTV = 124 L, RSE 11%
    lka_rtv    <- fixed(log(0.898))
    label("Ritonavir first-order absorption rate, ka (1/h, FIXED)")                  # Table 2: ka RTV (per hour) fixed = 0.898 (final estimate of separate RTV model in Table 1)
    lalag_rtv  <- log(1.05)
    label("Ritonavir absorption lag time, Tlag (h)")                                 # Table 2: Lag T RTV = 1.05 h, RSE 1%

    # ============================================================
    # Direct sigmoidal-Emax inhibition of ATV apparent clearance by
    # ritonavir plasma concentration:
    #   CL/F_ATV(t) = CL0/F_ATV * (1 - I(t))
    #   I(t)        = Imax * C_RTV / (IC50 + C_RTV)
    # Schipani 2013 Methods, "Population PK Modelling" section
    # (sigmoid function describing I(t)).
    # ============================================================
    imax  <- 0.988
    label("Maximum fractional inhibition of ATV CL/F by ritonavir (unitless)")       # Table 2: I max = 0.988, RSE 1%
    ic50  <- 0.221
    label("Ritonavir plasma concentration producing 50% of Imax on ATV CL/F (mg/L)") # Table 2: IC50 = 0.221 mg/L, RSE 13%

    # ============================================================
    # IIV - log-normal; omega^2 = log(1 + CV^2). Reported CV% from
    # Schipani 2013 Table 2 (final combined model).
    # The paper reports IIV only on V/F (both drugs) and on CL/F (RTV
    # only); IIV on CL/F ATV was tested and removed because it
    # destabilised the model (Schipani 2013 Discussion: "The addition
    # of IIV on ATV CL/F contributed to the instability of the model.
    # ... thus IIV on ATV CL/F was not included in the final model.").
    # ============================================================
    etalvc       ~ log(1 + 0.53^2)
    # Table 2: IIV V/F ATV = 53%  ->  omega^2 = log(1 + 0.53^2)

    # Correlated IIV block for ritonavir CL/F and V/F.
    # Covariance = rho * sqrt(var_cl * var_vc) where
    # var_cl = log(1 + 0.77^2), var_vc = log(1 + 0.73^2),
    # rho = 0.75 (Table 2: Correlation (CL_RTV, V_RTV) = 0.75).
    etalcl_rtv + etalvc_rtv ~ c(
      log(1 + 0.77^2),
      0.75 * sqrt(log(1 + 0.77^2) * log(1 + 0.73^2)),
      log(1 + 0.73^2)
    )
    # Table 2: IIV CL/F RTV = 77%, IIV V/F RTV = 73%, Correlation (CL_RTV, V_RTV) = 0.75

    # ============================================================
    # Residual error - proportional only in the final combined model
    # (Schipani 2013 Results, "Simultaneous Combined ATV-RTV Model":
    # "The residual variability was best described by a proportional
    # structure for ATV and RTV.").
    # ============================================================
    propSd      <- 0.63
    label("Atazanavir proportional residual error (fraction)")                       # Table 2: Proportional ATV = 63%, RSE 18%
    propSd_rtv  <- 0.73
    label("Ritonavir proportional residual error (fraction)")                        # Table 2: Proportional RTV = 73%, RSE 5%
  })

  model({
    # ------------------------------------------------------------
    # Individual ritonavir PK parameters. Computed first so the
    # ritonavir plasma concentration is available for the
    # atazanavir CL/F inhibition term below.
    # ------------------------------------------------------------
    cl_rtv   <- exp(lcl_rtv + etalcl_rtv)
    vc_rtv   <- exp(lvc_rtv + etalvc_rtv)
    ka_rtv   <- exp(lka_rtv)
    alag_rtv_h <- exp(lalag_rtv)

    # ------------------------------------------------------------
    # Ritonavir-driven sigmoidal-Emax inhibition of ATV CL/F.
    # In the absence of ritonavir (central_rtv = 0) the inhibition
    # evaluates to zero and ATV CL/F equals its unboosted typical
    # value exp(lcl) = 16.6 L/h.
    # ------------------------------------------------------------
    crtv     <- central_rtv / vc_rtv
    inhib    <- imax * crtv / (ic50 + crtv)

    # ------------------------------------------------------------
    # Individual atazanavir PK parameters. ATV CL/F is the only
    # structural parameter without an eta in the final combined
    # model (Schipani 2013 Discussion).
    # ------------------------------------------------------------
    cl       <- exp(lcl) * (1 - inhib)
    vc       <- exp(lvc + etalvc)
    ka       <- exp(lka)
    alag_atv_h <- exp(ltlag)

    # ------------------------------------------------------------
    # ODE system. Atazanavir occupies depot + central; ritonavir
    # occupies a parallel depot_rtv + central_rtv. The two drugs
    # are linked only via the inhibition term on ATV CL.
    # ------------------------------------------------------------
    d/dt(depot)        <- -ka * depot
    d/dt(central)      <-  ka * depot - (cl / vc) * central

    d/dt(depot_rtv)    <- -ka_rtv * depot_rtv
    d/dt(central_rtv)  <-  ka_rtv * depot_rtv - (cl_rtv / vc_rtv) * central_rtv

    # Absorption lag times (Tlag) for each depot.
    alag(depot)     <- alag_atv_h
    alag(depot_rtv) <- alag_rtv_h

    # ------------------------------------------------------------
    # Observation variables and proportional residual error.
    # Cc      = atazanavir plasma concentration (mg/L)
    # Cc_rtv  = ritonavir plasma concentration (mg/L)
    # ------------------------------------------------------------
    Cc     <- central     / vc
    Cc_rtv <- central_rtv / vc_rtv

    Cc     ~ prop(propSd)
    Cc_rtv ~ prop(propSd_rtv)
  })
}
