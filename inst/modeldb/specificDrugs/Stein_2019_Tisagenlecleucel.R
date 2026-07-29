Stein_2019_Tisagenlecleucel <- function() {
  description <- "Cellular kinetic model for tisagenlecleucel CAR-T cells in pediatric and young adult patients with relapsed or refractory B-cell acute lymphoblastic leukemia (Stein 2019). Single-infusion expansion-then-biexponential-decline analytical model: transgene levels grow exponentially at rate rho up to Tmax, after which effector cells decline at rate alpha and a fraction FB transitions to memory cells declining at rate beta."
  reference <- "Stein AM, Grupp SA, Levine JE, et al. Tisagenlecleucel Model-Based Cellular Kinetic Analysis of Chimeric Antigen Receptor-T Cells. CPT Pharmacometrics Syst Pharmacol. 2019;8(5):285-295. doi:10.1002/psp4.12388"
  vignette <- "Stein_2019_Tisagenlecleucel"
  units <- list(
    time          = "day",
    dosing        = "transgene copies/ug genomic DNA",
    concentration = "transgene copies/ug genomic DNA"
  )

  covariateData <- list(
    # No covariates were retained in the final model. Stein 2019 reports that
    # "bootstrapping the full covariate model showed no statistically
    # significant impact of any covariate on the maximal concentration of
    # tisagenlecleucel" (Results, Covariate analysis), so the final model
    # equals the base model and exposes no covariate columns.
  )

  population <- list(
    n_subjects     = 90,
    n_studies      = 2,
    age_range      = "3-25 years",
    age_median     = "12 years",
    weight_range   = "14-140 kg",
    weight_median  = "39 kg",
    sex_female_pct = 50,
    race_ethnicity = c(White = 77, Asian = 9, "Other/unknown" = 14),
    disease_state  = "Pediatric and young adult patients with relapsed or refractory B-cell acute lymphoblastic leukemia (r/r B-ALL).",
    dose_range     = "Median (range) 3.1e6 (0.2-5.4e6) CAR-positive viable T cells/kg in patients <=50 kg; total dose 1.0e8 (0.03-2.6e8) cells in patients >50 kg.",
    regions        = "Global; ELIANA (NCT02435849, 10 countries) and ENSIGN (NCT02228096, US).",
    notes          = "Stein 2019 Table 2 baseline demographics. Down syndrome 8%; previous stem cell transplant 57%; lymphodepleting chemotherapy with fludarabine 94%; received tocilizumab 36%; received corticosteroids 26%. Population summary refers to the 90 patients pooled in the model-based analysis."
  )

  ini({
    # Structural parameters from Stein 2019 Table 1 (final model = base model;
    # no covariates retained). Log-transformed because each is constrained
    # positive. Comedication effects Ftoci = 1.2 and Fster = 1.0 (Table 1)
    # are not implemented in this model file: see vignette Assumptions and
    # deviations for the rationale.
    lfoldx <- log(3900);   label("Fold expansion of tisagenlecleucel from baseline (foldx, unitless)")    # Stein 2019 Table 1
    ltmax  <- log(9.3);    label("Time to maximal expansion (Tmax, days)")                                # Stein 2019 Table 1
    lcmax  <- log(24000);  label("Maximum transgene-copy level after expansion (Cmax, copies/ug)")        # Stein 2019 Table 1
    lalpha <- log(0.16);   label("Rapid contraction rate constant (alpha, 1/day)")                        # Stein 2019 Table 1
    lfb    <- log(0.0079); label("Fraction of transgene copies surviving from peak into slow decline (FB, unitless)")  # Stein 2019 Table 1
    lbeta  <- log(0.0032); label("Slow terminal decline rate constant (beta, 1/day)")                     # Stein 2019 Table 1

    # Inter-individual variability. Stein 2019 Table 1 reports the random-effect
    # variances on the log-transformed parameters (Monolix omega^2); they are
    # entered here verbatim.
    etalfoldx ~ 2.4    # Stein 2019 Table 1 (Random effect foldx)
    etaltmax  ~ 0.38   # Stein 2019 Table 1 (Random effect Tmax)
    etalcmax  ~ 0.65   # Stein 2019 Table 1 (Random effect Cmax)
    etalalpha ~ 0.91   # Stein 2019 Table 1 (Random effect alpha)
    etalfb    ~ 0.8    # Stein 2019 Table 1 (Random effect FB)
    etalbeta  ~ 0.86   # Stein 2019 Table 1 (Random effect beta)

    # Residual error. Stein 2019 used a "constant error model in log-scale"
    # (additive on log[transgene]) which is equivalent to a proportional
    # residual error in linear space.
    propSd <- 0.56; label("Proportional residual error (transgene copies/ug, fraction)")  # Stein 2019 Table 1 (Residual error a)
  })

  model({
    # Individual parameters (log-normal IIV)
    foldx <- exp(lfoldx + etalfoldx)
    tmax  <- exp(ltmax  + etaltmax)
    cmax  <- exp(lcmax  + etalcmax)
    alpha <- exp(lalpha + etalalpha)
    FB    <- exp(lfb    + etalfb)
    beta  <- exp(lbeta  + etalbeta)

    # Derived rate constants (Stein 2019 Figure 2 Definitions)
    rho <- log(foldx) / tmax            # Expansion rate (1/day)
    R0  <- cmax / foldx                 # Baseline transcript level at t = 0 (copies/ug)
    AA  <- (1 - FB) * cmax              # Amplitude of fast decline at Tmax
    BB  <- FB * cmax                    # Amplitude of slow decline at Tmax

    # Analytical solution of the compartmental model (Stein 2019 Figure 2 and
    # Eq. 1; supplement section "Analytical solution for the compartmental
    # model"). Effector cells expand exponentially up to t = Tmax and then
    # decline biexponentially as the effector pool contracts (alpha) and a
    # fraction FB of cells transitions into the memory-like pool that decays
    # slowly (beta). The observed signal is the sum E + M.
    if (t < 0) {
      Cc <- R0
    } else if (t < tmax) {
      Cc <- R0 * exp(rho * t)
    } else {
      Cc <- AA * exp(-alpha * (t - tmax)) + BB * exp(-beta * (t - tmax))
    }

    Cc ~ prop(propSd)
  })
}
