Hayashi_2007_omalizumab <- function() {
  description <- "Mechanism-based binding population PK/PD model for omalizumab and IgE in Japanese atopic-asthma patients (Hayashi 2007). Three serum entities (free omalizumab, free IgE, and the omalizumab-IgE complex) each carry their own clearance and volume of distribution and are coupled through instantaneous-equilibrium binding (law of mass action) with a concentration-dependent dissociation constant. Body weight modifies omalizumab CL and Vd; baseline IgE modifies IgE CL and IgE production rate. Subcutaneous absorption is first-order. Disposition parameters are reported as apparent (divided by SC bioavailability f). Three observed quantities: total omalizumab (ug/mL), total IgE (ng/mL), and free IgE (ng/mL)."
  reference <- "Hayashi N, Tsukamoto Y, Sallas WM, Lowe PJ. A mechanism-based binding model for the population pharmacokinetics and pharmacodynamics of omalizumab. Br J Clin Pharmacol. 2007;63(5):548-561. doi:10.1111/j.1365-2125.2006.02803.x (PMID 17096680)."
  vignette <- "Hayashi_2007_omalizumab"
  units <- list(
    time          = "day",
    dosing        = "mg",
    concentration = "ug/mL (total omalizumab); ng/mL (free and total IgE)"
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power covariate on apparent CL of free omalizumab (CL_X/f, exponent 0.911) and on apparent V_d of omalizumab and IgE (V_X/f = V_E/f, exponent 0.658). Both effects centred at 61.1 kg (Hayashi 2007 Table 3 footnote *).",
      source_name        = "Body weight"
    ),
    IGE = list(
      description        = "Baseline serum total IgE concentration (pretreatment)",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power covariate on apparent CL of free IgE (CL_E/f, exponent -0.281) and on apparent IgE production rate (P_E/f, exponent 0.657). Both effects centred at 482.4 ng/mL (Hayashi 2007 Table 3 footnote dagger). Also used as the initial value for the total-IgE state total_target at t = 0: with no drug aboard, total_target(0) = (IGE / MW_IgE) * V_E. Source paper reports IU/mL alongside ng/mL with conversion 1 IU/mL = 2.42 ng/mL (Hayashi 2007 Methods, page 549).",
      source_name        = "IgE0 (baseline IgE concentration)"
    )
  )

  population <- list(
    n_subjects       = 202L,
    n_studies        = 2L,
    study_names      = c("1101 (single-dose, healthy atopic Japanese)",
                         "1305 (multiple-dose, Japanese seasonal allergic rhinitis)"),
    n_observations   = 3192L,
    n_omalizumab_obs = 1037L,
    n_total_ige_obs  = 1191L,
    n_free_ige_obs   = 964L,
    age_range        = "Adults; demographic table of Hayashi 2007 (Table 1) reports body weight and baseline IgE only",
    weight_range     = "Study 1101: 51-79 kg (mean 62.5, SD 6.4); Study 1305: 42-101 kg (mean 60.5, SD 10.2)",
    weight_median    = "Reference 61.1 kg (Hayashi 2007 Table 3 footnote *)",
    sex_female_pct   = NA_real_,
    race_ethnicity   = c(White = 0, Black = 0, Asian = 100, Other = 0),
    race_notes       = "All 202 model-building subjects were Japanese (Hayashi 2007 Table 1).",
    disease_state    = "Healthy Japanese atopic volunteers (study 1101) and Japanese patients with seasonal allergic rhinitis (study 1305)",
    dose_range       = "Study 1101: single SC 75, 150, 300, or 375 mg; Study 1305: multiple-dose SC 150-375 mg every 2 or 4 weeks per the dosing table (Hayashi 2007 Table 2) keyed by body weight and baseline IgE",
    regions          = "Japan",
    baseline_ige_range = "Study 1101: 204-2143 ng/mL (mean 811, SD 473); Study 1305: 53-1316 ng/mL (mean 373, SD 317)",
    notes            = "Model-building dataset from two Japanese clinical studies (1101 and 1305). The published external validation against 531 White patients from studies 007/008/009 is not part of the model-building cohort and is referenced only for predictive evaluation (Hayashi 2007 Table 1)."
  )

  ini({
    # Structural parameters - reference 61.1 kg body weight, 482.4 ng/mL baseline IgE.
    # All CL and V are reported as apparent (divided by SC bioavailability f).
    # Time unit converted from h (paper) to d for nlmixr2lib convention: x_per_day = x_per_h * 24.
    lka          <- log(0.0200 * 24);  label("Apparent SC absorption rate constant for omalizumab (ka, 1/d; paper: 0.0200 /h)")          # Hayashi 2007 Table 3
    lcl          <- log(0.00732 * 24); label("Apparent CL of free omalizumab at 61.1 kg (CL_X/f, L/d; paper: 7.32 mL/h)")                # Hayashi 2007 Table 3
    ldcl_complex <- log(0.00586 * 24); label("Apparent excess CL of complex over free omalizumab (Delta_CL_C/f, L/d; paper: 5.86 mL/h)") # Hayashi 2007 Table 3
    lcl_ige      <- log(0.071 * 24);   label("Apparent CL of free IgE at 482.4 ng/mL baseline (CL_E/f, L/d; paper: 71.0 mL/h)")          # Hayashi 2007 Table 3
    lvc          <- log(5.9);          label("Apparent V_d of omalizumab and IgE at 61.1 kg (V_X/f = V_E/f, L; paper: 5900 mL)")          # Hayashi 2007 Table 3 (V_E/f assumed equal to V_X/f, footnote double-dagger)
    lvc_complex  <- log(3.63);         label("Apparent V_d of omalizumab-IgE complex (V_C/f, L; paper: 3630 mL)")                          # Hayashi 2007 Table 3
    lp_ige       <- log(0.1595 * 24);  label("Apparent rate of IgE production at 482.4 ng/mL baseline (P_E/f, nmol/d; paper: 30.3 ug/h = 0.1595 nmol/h)") # Hayashi 2007 Table 3 footnote dagger (30.3 ug/h corresponds to 0.159 nmol/h)
    lkd0         <- log(1.07);         label("Equilibrium dissociation constant at central = total_target (Kd0, nM)")                                 # Hayashi 2007 Table 3

    # Covariate effects (power form; Hayashi 2007 page 555 equations).
    e_wt_cl      <-  0.911;  label("Power exponent of body weight on CL_X/f (unitless)")     # Hayashi 2007 Table 3
    e_wt_vc      <-  0.658;  label("Power exponent of body weight on V_X/f (unitless)")      # Hayashi 2007 Table 3
    e_ige_cl_ige <- -0.281;  label("Power exponent of baseline IgE on CL_E/f (unitless)")    # Hayashi 2007 Table 3
    e_ige_p_ige  <-  0.657;  label("Power exponent of baseline IgE on P_E/f (unitless)")     # Hayashi 2007 Table 3
    alpha        <-  0.157;  label("Concentration-dependence exponent on Kd (unitless)")     # Hayashi 2007 Table 3

    # IIV - log-normal; omega^2 = log(CV^2 + 1) for log-normal parameters.
    # ka     CV 39.9% -> log(1 + 0.399^2) = 0.14773
    # CL_X   CV 20.3% -> log(1 + 0.203^2) = 0.04042
    # dCL_C  CV 34.9% -> log(1 + 0.349^2) = 0.11488
    # V_X    CV 13.0% -> log(1 + 0.130^2) = 0.01679
    # V_C    CV 25.0% -> log(1 + 0.250^2) = 0.06062
    # CL_E   CV 25.3% -> log(1 + 0.253^2) = 0.06205
    # P_E    CV 23.1% -> log(1 + 0.231^2) = 0.05197
    # corr(CL_E, P_E) = 0.968 -> cov = 0.968 * sqrt(0.06205 * 0.05197) = 0.05496
    etalka          ~ 0.14773                                                                            # Hayashi 2007 Table 3 (ka CV 39.9%)
    etalcl          ~ 0.04042                                                                            # Hayashi 2007 Table 3 (CL_X CV 20.3%)
    etaldcl_complex ~ 0.11488                                                                            # Hayashi 2007 Table 3 (Delta_CL_C CV 34.9%)
    etalvc          ~ 0.01679                                                                            # Hayashi 2007 Table 3 (V_X CV 13.0%; V_E/f = V_X/f shares this eta)
    etalvc_complex  ~ 0.06062                                                                            # Hayashi 2007 Table 3 (V_C CV 25.0%)
    etalcl_ige + etalp_ige ~ c(0.06205,
                          0.05496, 0.05197)                                                              # Hayashi 2007 Table 3 (CL_E CV 25.3%, P_E CV 23.1%, correlation 0.968)

    # Residual error - paper uses Y = F * exp(eps), eps ~ N(0, sigma^2) (Hayashi 2007 page 553).
    # For small sigma this is well approximated by a proportional model
    # Y ~ F + F * propSd * eps with propSd = sigma; encoded as prop() per nlmixr2 convention.
    propSd        <- 0.167; label("Proportional residual error on total omalizumab (fraction)")        # Hayashi 2007 Table 3 (omalizumab intra-individual CV 16.7%)
    propSd_totalIgE <- 0.211; label("Proportional residual error on total IgE (fraction)")               # Hayashi 2007 Table 3 (total IgE intra-individual CV 21.1%)
    propSd_freeIgE  <- 0.218; label("Proportional residual error on free IgE (fraction)")                # Hayashi 2007 Table 3 (free IgE intra-individual CV 21.8%)
  })

  model({
    # ------------------------------------------------------------------
    # Constants - molecular weights from Hayashi 2007 Methods page 552.
    # 1 nM * MW_kDa = MW_kDa ng/mL (since 1 nmol/L * 1 g/mol = 1e-9 g/L
    # = 1e-6 mg/L = 1 ng/mL when MW is taken in g/mol). For dose-side
    # mass-to-mole conversion, 1 mg = 1000 / MW_kDa nmol.
    # ------------------------------------------------------------------
    MWX <- 150  # omalizumab molecular weight (kDa = ng/nmol = g/mmol)
    MWE <- 190  # IgE molecular weight (kDa)

    # ------------------------------------------------------------------
    # 1. Individual parameters.
    # ------------------------------------------------------------------
    ka          <- exp(lka          + etalka)
    cl          <- exp(lcl          + etalcl)          * (WT  / 61.1)^e_wt_cl
    dcl_complex <- exp(ldcl_complex + etaldcl_complex)
    cl_complex  <- cl + dcl_complex
    cl_ige      <- exp(lcl_ige      + etalcl_ige)      * (IGE / 482.4)^e_ige_cl_ige
    p_ige       <- exp(lp_ige       + etalp_ige)       * (IGE / 482.4)^e_ige_p_ige
    vc          <- exp(lvc          + etalvc)          * (WT  / 61.1)^e_wt_vc
    v_ige       <- vc                                # V_E/f = V_X/f (Hayashi 2007 Table 3 footnote double-dagger)
    vc_complex  <- exp(lvc_complex  + etalvc_complex)
    kd0         <- exp(lkd0)

    # ------------------------------------------------------------------
    # 2. Concentration-dependent dissociation constant Kd
    #    (Hayashi 2007 "The dissociation constant" page 552):
    #    Kd = Kd0 * (central / total_target)^alpha. At t = 0, central = 0 so Kd = 0
    #    (0^alpha = 0 for alpha > 0) and no complex forms.
    # ------------------------------------------------------------------
    kd <- kd0 * (central / total_target)^alpha

    # ------------------------------------------------------------------
    # 3. Equilibrium-binding solution for the complex amount X_C.
    #    Mass-balance + law of mass action with K_d = C_fX * C_fE / C_C
    #    yields the quadratic
    #      X_C^2 - X_C * (central + total_target + Kd * V_X * V_E / V_C) + central * total_target = 0
    #    whose negative root (Hayashi 2007 page 552) is
    #      X_C = (1/2) * { S - sqrt(S^2 - 4 * central * total_target) }
    #    with S = central + total_target + Kd * V_X * V_E / V_C. The discriminant is
    #    >= 0 for all physical central, total_target, Kd, V's.
    # ------------------------------------------------------------------
    S   <- central + total_target + kd * vc * v_ige / vc_complex
    X_C <- 0.5 * (S - sqrt(S * S - 4 * central * total_target))

    # Free / complex concentrations (nM = nmol/L).
    C_fX <- (central - X_C) / vc
    C_fE <- (total_target - X_C) / v_ige
    C_C  <-  X_C       / vc_complex

    # ------------------------------------------------------------------
    # 4. ODE system (Hayashi 2007 "PK/PD structural model" page 552).
    #    State units: depot (mg), central, total_target, X_C (nmol).
    #    Rate units: ka * depot (mg/d), p_ige (nmol/d), cl * C (L/d * nM = nmol/d).
    #    The depot-to-central transfer converts mg to nmol via 1000 / MWX.
    # ------------------------------------------------------------------
    d/dt(depot) <- -ka * depot
    d/dt(central)  <-  ka * depot * (1000 / MWX) - cl     * C_fX - cl_complex * C_C
    d/dt(total_target)  <-  p_ige                - cl_ige * C_fE - cl_complex * C_C

    # Initial total IgE amount at t = 0 (no drug, so total IgE = free IgE):
    #   total_target(0) = (IGE_ng_per_mL / MWE_kDa) * V_E_L = nM * L = nmol.
    total_target(0) <- (IGE / MWE) * v_ige

    # ------------------------------------------------------------------
    # 5. Observation outputs in assay units.
    #    Total omalizumab in ug/mL = (C_fX + C_C) [nM] * MWX [ng/nmol] / 1000.
    #    Total IgE in ng/mL        = (C_fE + C_C) [nM] * MWE [ng/nmol].
    #    Free IgE in ng/mL         = C_fE [nM] * MWE [ng/nmol].
    # ------------------------------------------------------------------
    Cc       <- (C_fX + C_C) * MWX / 1000
    totalIgE <- (C_fE + C_C) * MWE
    freeIgE  <-  C_fE        * MWE

    Cc       ~ prop(propSd)
    totalIgE ~ prop(propSd_totalIgE)
    freeIgE  ~ prop(propSd_freeIgE)
  })
}
