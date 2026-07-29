Csajka_2005_ephedrine_caffeine <- function() {
  description <- "Mechanistic simultaneous population PK model for co-administered ephedrine, its N-demethylation metabolite norephedrine, and caffeine in healthy adults after single oral doses (Csajka 2005). Caffeine is described by a 1-compartment first-order-absorption model with a fractional decrease in apparent clearance during oral contraceptive therapy. Ephedrine uses a 1-compartment depot + central + cumulative-urine model with an absorption lag time, renal clearance, and saturable Michaelis-Menten conversion to norephedrine; norephedrine is carried as a pseudo-concentration state because its volume of distribution V_NE is unidentifiable, so the reported parameter is the compound Vmax/V_NE and the norephedrine elimination is first order. The interaction term reproduces the paper's indirect-action absorption model (equation 10b/10e final form): the caffeine amount in the absorption compartment depresses ephedrine ka by an asymptotic fraction d, with caffeine acting as the f(C) inhibitor on its own absorption-compartment amount. Parameter values are the pharmaceutical-formulation defaults from Table 3; herbal-formulation alternatives (bioavailability F_E,herbal = 0.78 instead of F_E,pharm = 0.59, plus a 22.2-min caffeine absorption lag) are documented in inline comments and can be applied by overriding lfdepot and ltlag_caf at simulation time."
  reference   <- "Csajka C, Haller CA, Benowitz NL, Verotta D. Mechanistic pharmacokinetic modelling of ephedrine, norephedrine and caffeine in healthy subjects. Br J Clin Pharmacol. 2005;59(3):335-345. doi:10.1111/j.1365-2125.2005.02254.x"
  vignette    <- "Csajka_2005_ephedrine_caffeine"

  units <- list(time = "min", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CONMED_BIRTHCONTROL = list(
      description        = "Oral hormonal contraceptive use indicator (1 = currently taking, 0 = not)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no oral contraceptive)",
      notes              = "Csajka 2005 Methods 'Caffeine pharmacokinetics' parameterises the OC effect as CL_C(OC) = CL_C * (1 - d_OC_CL * OC), with d_OC_CL = 0.54 (Table 3), so OC users have caffeine apparent clearance reduced by 54% (population mean 0.083 -> 0.038 L/min). Six of 24 subjects (25%; two from Study 1 and four from Study 2) were on oral contraceptives. The covariate affects only caffeine CL; no effect on ephedrine or norephedrine.",
      source_name        = "OC"
    )
  )

  covariatesDataExcluded <- list(
    URINE_PH = list(
      description = "Voiding-interval urine pH",
      units       = "pH units",
      type        = "continuous",
      notes       = "Csajka 2005 reports an inverse linear association between individual measurements of urine pH and empirical-Bayes individual estimates of ephedrine renal clearance (CL_R = 0.4723 - 0.0172 * pH, P = 0.013, Results 'Ephedrine and norephedrine pharmacokinetics'). This was an exploratory post-hoc regression on EBE estimates and was NOT incorporated as a covariate in the population model; the final model carries CL_RE as a fixed-effect population mean with no pH covariate. Documented here so the post-hoc finding is not lost to future users."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 24L,
    n_studies      = 2L,
    age_range      = "22-39 years (Study 1 25-38; Study 2 22-39)",
    weight_range   = "52-91.5 kg (Study 1 52-88.9; Study 2 58.6-91.5)",
    height_range   = "144-188 cm (Study 2 only; not reported for Study 1)",
    sex_female_pct = 58.3,
    race_ethnicity = c(White = 62.5, Black_or_African_American = 6.25, Asian_or_Pacific_Islander = 18.75, Hispanic_or_Latino = 12.5),
    disease_state  = "Healthy adults. 6 of 24 subjects (25%) were on concurrent oral contraceptive therapy. Race / ethnicity reported only for Study 2 (n=16: 10 Caucasian, 1 African American, 3 Asian / Pacific Islander, 2 Latino); Study 1 (n=8) race not reported.",
    dose_range     = "Study 1 (n=8): single oral dose of a commercial herbal supplement (Metabolift, two capsules) containing 17.3 mg ephedrine, 0.2 mg norephedrine, 5.3 mg pseudoephedrine, 0.42 mg norepseudoephedrine, and 175 mg caffeine. Study 2 (n=16): single oral 25 mg ephedrine sulphate (West-ward Pharmaceutical Corp.) or 200 mg caffeine sulphate alone and together; cross-over.",
    regions        = "USA (University of California, San Francisco)",
    notes          = "Csajka 2005 Table 1 demographics. 379 ephedrine, 352 norephedrine, 417 caffeine plasma samples and 40 urinary ephedrine collections (0-14 h Study 1; 0-24 h Study 2). NONMEM FOCE-INTERACTION with ADVAN6. Pharmaceutical-formulation parameters are the defaults in this model; the herbal formulation differs in caffeine absorption lag (22.2 min vs 0) and apparent ephedrine bioavailability (F_E,herbal 0.78 vs F_E,pharm 0.59) per Table 3."
  )

  ini({
    # =========================================================
    # Caffeine 1-compartment PK with first-order absorption.
    # Csajka 2005 Table 3, caffeine section. Apparent parameters
    # CL_C / F_C and V_C / F_C: oral bioavailability F_C is
    # implicit in the F-adjusted CL and V. Time in minutes,
    # amounts in mg, volumes in L; the table is reported in
    # L / min and L; concentration units mg / L convert directly
    # to ug / mL (matches the paper's reported concentration
    # range of 0 to 8470 ug / L for plasma caffeine).
    # =========================================================
    lka_caf <- log(0.064);  label("Caffeine first-order absorption rate ka_C (1/min)")           # Table 3 (ka_C = 0.064 1/min)
    lcl_caf <- log(0.083);  label("Caffeine apparent oral clearance CL_C/F_C (L/min)")           # Table 3 (CL_C/F_C = 0.083 L/min, no OC)
    lvc_caf <- log(38.6);   label("Caffeine apparent central volume V_C/F_C (L)")                # Table 3 (V_C/F_C = 38.6 L)

    # Fractional decrease in caffeine apparent CL during oral
    # contraceptive therapy. Csajka 2005 Table 3 (d_OC_CL =
    # 0.54). Multiplicative form: CL_C(OC) = CL_C * (1 - d_OC).
    e_conmed_birthcontrol_cl_caf <- 0.54;  label("Fractional decrease in caffeine CL_C/F_C during oral contraceptive therapy (unitless)") # Table 3 (d_OC_CL = 0.54)

    # =========================================================
    # Ephedrine 1-compartment PK with first-order absorption,
    # lag time, apparent renal clearance, and parallel saturable
    # metabolism to norephedrine. Csajka 2005 Table 3, ephedrine
    # section. The paper labels CL_RE as renal clearance because
    # ephedrine is eliminated primarily by renal excretion; the
    # MM metabolic pathway to norephedrine is parameterised
    # separately via lvmax_neph / lkm. F_E is implicit in the
    # apparent CL and V, so CL_RE / F_E and V_E / F_E are the
    # reported values; lfdepot is a separate explicit
    # bioavailability factor for distinguishing the pharmaceutical
    # (default) and herbal formulations.
    # =========================================================
    lka         <- log(0.036);  label("Ephedrine first-order absorption rate ka_E without caffeine (1/min)") # Table 3 (ka_E = 0.036 1/min, no caffeine)
    lcl_renal   <- log(0.34);   label("Ephedrine apparent renal clearance CL_RE/F_E (L/min)")               # Table 3 (CL_RE/F_E = 0.34 L/min)
    lvc         <- log(181);    label("Ephedrine apparent central volume V_E/F_E (L)")                      # Table 3 (V_E/F_E = 181 L)
    ltlag       <- log(16.7);   label("Ephedrine absorption lag time (min)")                                # Table 3 (lag time = 16.7 min)

    # Ephedrine bioavailability anchor F_E. Default = F_E,pharm =
    # 0.59 (pharmaceutical formulation). For herbal-formulation
    # simulation override with lfdepot <- log(0.78) per Table 3
    # (F_E,herbal). The IIV on lfdepot is set to fixed(0) because
    # the paper reports "No significant intersubject variability
    # was detected for these parameters and the variance was fixed
    # to zero" (Results, Influence of formulation).
    lfdepot     <- fixed(log(0.59)); label("Ephedrine apparent bioavailability F_E,pharm (unitless)")        # Table 3 (F_E,pharm = 0.59)

    # =========================================================
    # Norephedrine metabolite parameters. The MM rate of
    # ephedrine -> norephedrine and norephedrine elimination.
    # Vmax/V_NE is a compound parameter because the norephedrine
    # volume V_NE cannot be identified (norephedrine is not
    # administered). The state central_neph therefore carries the
    # norephedrine pseudo-concentration directly (= A5 / V_NE) so
    # that (Vmax/V_NE) * MM_factor has the proper concentration-
    # rate units. See model() block for the explicit ODE form.
    # =========================================================
    lvmax_neph  <- log(1.96e-4); label("Norephedrine compound Vmax/V_NE for ephedrine -> norephedrine (mg / (min L))")  # Table 3 (Vmax/V_NE = 1.96e-4 mg / (min L))
    lkm         <- log(2.77);    label("Michaelis-Menten Km for ephedrine -> norephedrine, on ephedrine central amount (mg)") # Table 3 (Km = 2.77 mg)
    lkel_neph   <- log(0.037);   label("Norephedrine first-order elimination rate ke_NE (1/min)")  # Table 3 (ke_NE = 0.037 1/min) -- table heading is mislabelled as 'L/h'; the prose text states '0.037 1 min-1', which is the value used here

    # =========================================================
    # Caffeine effect on ephedrine absorption (Csajka 2005
    # equation 10b / 10e, the final indirect-action model).
    # ka_E is depressed when caffeine is present in its
    # absorption compartment: ka_E(t) = ka_E * (1 - d_caf *
    # depot_caf / (a50_caf + depot_caf)), where d_caf is the
    # asymptotic decrease (d = exp(theta_e)/(1 + exp(theta_e)))
    # and a50_caf is the amount of caffeine in the absorption
    # compartment producing a 50% fractional decrease. Reported
    # in Table 3 as theta_e = 19.5 (giving d = 0.99) and
    # ka_E,50 = 31.4 mg. No IIV reported for these parameters.
    # =========================================================
    la50_caf    <- log(31.4); label("Caffeine amount in absorption compartment at half-max ephedrine ka_E decrease (mg)")  # Table 3 (ka_E,50 = 31.4 mg)
    d_caf       <- fixed(0.99); label("Asymptotic fractional decrease in ephedrine ka_E due to caffeine in absorption compartment (unitless, 0-1)")  # Table 3 (d = 0.99 = expit(theta_e = 19.5))

    # =========================================================
    # Inter-individual variability. Csajka 2005 Table 3 reports
    # %CV; converted via omega^2 = log(CV^2 + 1) (log-normal
    # exponential random-effect on typical value). Only
    # parameters for which the paper reports an IIV estimate are
    # included.
    # =========================================================
    etalka_caf  ~ 0.2231       # Table 3 caffeine ka_C %CV = 50; omega^2 = log(0.50^2 + 1) = 0.2231
    etalcl_caf  ~ 0.1346       # Table 3 caffeine CL_C %CV = 38; omega^2 = log(0.38^2 + 1) = 0.1346
    etalvc_caf  ~ 0.03922      # Table 3 caffeine V_C  %CV = 20; omega^2 = log(0.20^2 + 1) = 0.03922

    etalka      ~ 0.3266       # Table 3 ephedrine ka_E %CV = 63; omega^2 = log(0.63^2 + 1) = 0.3266
    etalcl_renal ~ 0.01203     # Table 3 ephedrine CL_RE %CV = 11; omega^2 = log(0.11^2 + 1) = 0.01203
    etalvc      ~ 0.03568      # Table 3 ephedrine V_E  %CV = 19; omega^2 = log(0.19^2 + 1) = 0.03568
    etaltlag    ~ 0.1983       # Table 3 ephedrine lag %CV = 47; omega^2 = log(0.47^2 + 1) = 0.1983

    etalvmax_neph ~ 0.01203    # Table 3 norephedrine Vmax/V_NE %CV = 11; omega^2 = log(0.11^2 + 1) = 0.01203
    etalkm        ~ 0.1689     # Table 3 norephedrine Km %CV = 43; omega^2 = log(0.43^2 + 1) = 0.1689
    etalkel_neph  ~ 0.1840     # Table 3 norephedrine ke_NE %CV = 45; omega^2 = log(0.45^2 + 1) = 0.1840

    # =========================================================
    # Residual error. Csajka 2005 Table 3 reports proportional
    # residual error as %CV; SD is %CV / 100 on the proportional
    # scale.
    # =========================================================
    propSd      <- 0.17;  label("Ephedrine plasma proportional residual SD (fraction)")     # Table 3 (ephedrine plasma s = 17%)
    propSd_caf  <- 0.17;  label("Caffeine plasma proportional residual SD (fraction)")      # Table 3 (caffeine plasma s = 17%)
    propSd_neph <- 0.15;  label("Norephedrine plasma proportional residual SD (fraction)")  # Table 3 (norephedrine plasma s = 15%)
  })

  model({
    # ------------------------------------------------------------
    # 1. Individual PK parameters. Exponential random-effect on
    # typical value.
    # ------------------------------------------------------------
    ka_caf  <- exp(lka_caf  + etalka_caf)
    cl_caf  <- exp(lcl_caf  + etalcl_caf) * (1 - e_conmed_birthcontrol_cl_caf * CONMED_BIRTHCONTROL)
    vc_caf  <- exp(lvc_caf  + etalvc_caf)

    ka      <- exp(lka      + etalka)
    cl_renal <- exp(lcl_renal + etalcl_renal)
    vc      <- exp(lvc      + etalvc)
    tlag    <- exp(ltlag    + etaltlag)

    vmax_neph <- exp(lvmax_neph + etalvmax_neph)
    km        <- exp(lkm        + etalkm)
    kel_neph  <- exp(lkel_neph  + etalkel_neph)

    a50_caf <- exp(la50_caf)

    # ------------------------------------------------------------
    # 2. Caffeine indirect-action effect on ephedrine absorption
    # rate. ka_E(t) is depressed by an asymptotic fraction d_caf
    # when caffeine is present in its absorption compartment
    # (Csajka 2005 equation 10b / 10e, final model). At
    # depot_caf >> a50_caf the absorption is suppressed to (1 -
    # d_caf) * ka, and at depot_caf -> 0 the absorption recovers
    # to the unperturbed ka.
    # ------------------------------------------------------------
    ka_eff <- ka * (1 - d_caf * depot_caf / (a50_caf + depot_caf))

    # ------------------------------------------------------------
    # 3. Caffeine ODE block. 1-compartment first-order absorption
    # + first-order elimination. depot_caf and central_caf hold
    # amounts in mg.
    # ------------------------------------------------------------
    d/dt(depot_caf)   <- -ka_caf * depot_caf
    d/dt(central_caf) <-  ka_caf * depot_caf - (cl_caf / vc_caf) * central_caf

    # ------------------------------------------------------------
    # 4. Ephedrine ODE block. 1-compartment first-order absorption
    # with caffeine-modulated ka and renal elimination to a
    # cumulative-urine compartment. Lag time and bioavailability
    # are applied via alag(depot) and f(depot).
    #
    # The norephedrine metabolic pathway is a parallel-formation
    # readout only -- no metabolic depletion term is subtracted
    # from central. The paper's discussion (Csajka 2005, Discussion
    # paragraph 3) states the conversion to norephedrine is "minor
    # compared with the renal elimination of ephedrine and, even
    # after repeated drug intake, significant drug accumulation
    # would not be expected"; equivalently the published Vmax/V_NE
    # determines the formation rate of norephedrine but the
    # absolute Vmax that would deplete ephedrine is not identified
    # from the data (V_NE is unknown). The standard mechanistic
    # simplification is to treat the metabolite formation as a
    # parallel readout of the central state rather than as an
    # additional depletion arm.
    # ------------------------------------------------------------
    mm_factor <- central / (km + central)
    d/dt(depot)   <- -ka_eff * depot
    d/dt(central) <-  ka_eff * depot - cl_renal * (central / vc)
    d/dt(urine)   <-  cl_renal * (central / vc)

    # ------------------------------------------------------------
    # 5. Norephedrine pseudo-concentration ODE. central_neph =
    # A5 / V_NE carries the concentration directly (V_NE is
    # unidentifiable from the data). The formation-rate factor
    # (Vmax/V_NE) * MM_factor has units of mg / (min L) = the
    # rate of change of the norephedrine concentration.
    # ------------------------------------------------------------
    d/dt(central_neph) <- vmax_neph * mm_factor - kel_neph * central_neph

    # ------------------------------------------------------------
    # 6. Lag time and bioavailability for the ephedrine depot.
    # The pharmaceutical formulation (default) carries F_E,pharm
    # = 0.59; herbal users should override lfdepot.
    # ------------------------------------------------------------
    alag(depot) <- tlag
    f(depot)    <- exp(lfdepot)

    # ------------------------------------------------------------
    # 7. Observations.
    # Cc       -- ephedrine plasma concentration (mg/L; multiply
    #             by 1000 to convert to ug/L = ng/mL).
    # Cc_caf   -- caffeine plasma concentration (mg/L).
    # Cc_neph  -- norephedrine plasma pseudo-concentration (mg/L).
    # The cumulative urinary ephedrine amount is available as the
    # state `urine` (mg) for downstream comparison against 0-14h
    # / 0-24h urine totals; no residual error is applied here
    # since the model file is structural-only.
    # ------------------------------------------------------------
    Cc      <- central     / vc
    Cc_caf  <- central_caf / vc_caf
    Cc_neph <- central_neph

    Cc      ~ prop(propSd)
    Cc_caf  ~ prop(propSd_caf)
    Cc_neph ~ prop(propSd_neph)
  })
}
