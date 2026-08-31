Proctor_2026_durvalumab_poppk <- function() {
  description <- "Empirical two-compartment popPK model of durvalumab (anti-PD-L1) with parallel linear and Michaelis-Menten clearance and sigmoidal time-dependent decline in linear clearance. This is the comparator model that Proctor 2026 rebuilt in RxODE and simulated against its own two-pore PBPK model (Figures 4 and 5); the parameter values originate in the Baverel 2018 durvalumab popPK analysis and are transcribed from Proctor 2026 Table S3."
  reference <- "Proctor JR, Wong H. Albumin Levels Are Predictive of Cachexia-Induced Time-Dependent Clearance of Therapeutic Antibodies: A Physiologically Based Pharmacokinetic Model of Durvalumab. CPT Pharmacometrics Syst Pharmacol. 2026;15(1):e70185. doi:10.1002/psp4.70185 (Table S3, row 'Durvalumab'; parameters originally reported by Baverel PG, Dubois VFS, Jin CY, et al. Population Pharmacokinetics of Durvalumab in Cancer Patients and Association With Longitudinal Biomarkers of Disease Status. Clin Pharmacol Ther. 2018;103(4):631-642. doi:10.1002/cpt.982)"
  vignette <- "Proctor_2026_durvalumab_cachexia"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. verified = TRUE against Proctor 2026 Table S3, whose
  # footnote defines Vc as the central and Vp as the peripheral volume.
  compartmentData <- list(
    central     = list(analyte = "durvalumab", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "durvalumab", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = NA_integer_,
    n_studies      = NA_integer_,
    disease_state  = "Adults with advanced solid tumours receiving durvalumab monotherapy.",
    dose_range     = "Simulated at 10 mg/kg IV every 2 weeks for 60 weeks in Proctor 2026 Methods 2.2.",
    notes          = paste(
      "PROVENANCE: this is a SECONDARY transcription. Proctor 2026 Table S3 lists",
      "the structural parameters for ten immune-checkpoint inhibitors taken from",
      "their original popPK publications; the durvalumab row is from Baverel 2018,",
      "which is not on disk for this extraction. Proctor 2026 states that it",
      "rebuilt this model in RxODE and simulated it as the empirical comparator",
      "for its PBPK predictions, so the model as packaged here is exactly the one",
      "Proctor 2026 ran. Re-extract directly from Baverel 2018 when that paper is",
      "acquired: Table S3 carries only the nine structural parameters below and",
      "none of Baverel 2018's covariate model, IIV, residual error, or the",
      "non-parametric confidence intervals that Proctor 2026 bootstrapped in",
      "Figure 5."
    ),
    scope_note     = paste(
      "Typical-value model only: Table S3 reports no IIV and no residual error, so",
      "both are fixed at zero. The other nine checkpoint inhibitors in Table S3",
      "are NOT packaged as model files -- they are other authors' models",
      "transcribed for a figure, and three of them already have primary-source",
      "extractions in this library (Lindauer_2017_pembrolizumab,",
      "Bajaj_2017_nivolumab_ddmore, Hwang_2022_tremelimumab). Their Table S3",
      "time-dependent-clearance parameters are reproduced in the vignette for the",
      "Figure 6 replication only."
    )
  )

  ini({
    # Structural parameters, Proctor 2026 Table S3 row "Durvalumab". All fixed:
    # they are transcribed point estimates, not re-estimated here.
    lcl   <- fixed(log(0.249))  ; label("Baseline linear clearance CL0 (L/day)")            # Table S3 CL0 0.249 L/day
    lvc   <- fixed(log(3.51))   ; label("Central volume of distribution (L)")               # Table S3 Vc 3.51 L
    lq    <- fixed(log(0.477))  ; label("Inter-compartmental clearance (L/day)")            # Table S3 Q 0.477 L/day
    lvp   <- fixed(log(3.56))   ; label("Peripheral volume of distribution (L)")            # Table S3 Vp 3.56 L
    lvmax <- fixed(log(0.744))  ; label("Maximal rate of non-linear clearance (mg/day)")    # Table S3 Vmax 0.744 mg/day
    lkm   <- fixed(log(0.452))  ; label("Concentration at half-maximal non-linear clearance (mg/L)")  # Table S3 KM 0.452 mg/L

    # Sigmoidal time-dependent decline in linear clearance, Eq 2:
    #   CL(t) = CL0 * exp(Imax * t^gamma / (TI50^gamma + t^gamma))
    # Imax is the LOGARITHM of the maximal fold-change in clearance (paper
    # Methods 2.2), so it is negative for a declining clearance and is NOT
    # log-transformed here. This form is confirmed numerically: at t = 420 days
    # it gives exp(-0.185 * 0.8336) - 1 = -14.3%, matching the paper's statement
    # that the popPK model predicts a 14% decline in CL after 60 weeks.
    cl_time_max    <- fixed(-0.185)     ; label("Log of the maximal fold-change in linear clearance over time (Imax)")  # Table S3 Imax -0.185
    lcl_t50   <- fixed(log(173.1)) ; label("Time of half-maximal change in clearance TI50 (day)")                  # Table S3 TI50 173.1 day
    lcl_time_hill <- fixed(log(1.817)) ; label("Hill slope of the time-dependent clearance function (gamma)")           # Table S3 gamma 1.817

    # Table S3 reports no residual error; fixed at zero per the standing policy
    # for unreported RUV (documented in the vignette Errata).
    propSd <- fixed(0)          ; label("Proportional residual error (fraction)")
  })

  model({
    cl0            <- exp(lcl)
    vc             <- exp(lvc)
    q              <- exp(lq)
    vp             <- exp(lvp)
    vmax           <- exp(lvmax)
    km             <- exp(lkm)
    cl_t50    <- exp(lcl_t50)
    cl_time_hill  <- exp(lcl_time_hill)

    # Time-dependent linear clearance (Eq 2). t is in days, matching TI50.
    # t is clamped at zero because a negative base raised to the non-integer
    # Hill exponent is undefined; clamping returns CL0 before the first dose,
    # which is the correct pre-treatment value rather than a silent fudge.
    tdcl <- max(t, 0)
    cl <- cl0 * exp(cl_time_max * tdcl^cl_time_hill /
                      (cl_t50^cl_time_hill + tdcl^cl_time_hill))

    k12 <- q / vc
    k21 <- q / vp

    Cc <- central / vc

    # Parallel linear and Michaelis-Menten elimination from the central
    # compartment; Vmax is already an amount rate (mg/day) so it is not scaled.
    d/dt(central)     <- -cl / vc * central - k12 * central + k21 * peripheral1 -
      vmax * Cc / (km + Cc)
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc ~ prop(propSd)
  })
}
