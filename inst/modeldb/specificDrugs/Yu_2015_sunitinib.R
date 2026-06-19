Yu_2015_sunitinib <- function() {
  description <- paste0(
    "Integrated semi-physiological population PK model for oral sunitinib ",
    "and its equipotent active metabolite N-desethyl sunitinib (SU12662) in ",
    "adult cancer patients (n = 70 across three studies). Sunitinib is ",
    "absorbed first-order into a hypothetical hepatic enzyme compartment ",
    "that sits algebraically in equilibrium with the sunitinib central ",
    "compartment via hepatic blood flow Qh (fixed at 80 L/h for a 70 kg ",
    "subject). Clearance CL of sunitinib acts at the enzyme site (Cliv = ",
    "(ka * depot + Qh / Vc * central) / (Qh + CL)); fraction fm = 0.21 ",
    "(fixed from Houk 2009) of the cleared sunitinib appears as SU12662 ",
    "input into the metabolite central compartment, the rest is true ",
    "(non-SU12662) sunitinib clearance. SU12662 follows a 2-compartment ",
    "disposition with its own central and peripheral volumes and inter-",
    "compartmental clearance. Body-weight allometric scaling with fixed ",
    "exponents 0.75 (clearance / flow: CL_sun, Qh, CL_SU12662, Qi_SU12662) ",
    "and 1.0 (volumes: Vc_sun, Vc_SU12662, Vp_SU12662) is applied a priori ",
    "with reference WT = 70 kg. IIV is modelled with a 4 x 4 OMEGA BLOCK on ",
    "Vc_sun, Vc_SU12662, CL_SU12662, CL_sun with the paper's reported ",
    "correlations 0.48, 0.45, 0.53 and remaining off-diagonals fixed to ",
    "zero. Per-study residual proportional error from Yu 2015 Table 2 is ",
    "simplified to a single propSd / propSd_su12662 pair populated from the ",
    "Study 1 sigma^2 estimates (the largest cohort: 50 patients, 703 ",
    "samples); the smaller Study 2 + 3 residuals are noted in the vignette."
  )
  reference <- paste(
    "Yu H, Steeghs N, Kloth JSL, de Wit D, van Hasselt JGC, van Erp NP,",
    "Beijnen JH, Schellens JHM, Mathijssen RHJ, Huitema ADR.",
    "Integrated semi-physiological pharmacokinetic model for both sunitinib",
    "and its active metabolite SU12662.",
    "Br J Clin Pharmacol. 2015;79(5):809-819.",
    "doi:10.1111/bcp.12550.",
    sep = " "
  )
  vignette <- "Yu_2015_sunitinib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at baseline; used for fixed-exponent allometric scaling with reference 70 kg. Yu 2015 Methods: exponent 0.75 on CL / Qh / Qi and 1.0 on Vc / Vp (Eqs 3-4). About 6% of subjects in the paper had no WT and were imputed to the population mean of 70 kg; the model assumes a finite positive WT is supplied per subject.",
      source_name        = "WT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 70L,
    n_studies      = 3L,
    age_range      = "adult cancer patients (specific range not extractable from main text)",
    weight_range   = "39-157 kg",
    weight_median  = "82 kg",
    sex_female_pct = NA_real_,
    race_ethnicity = NULL,
    disease_state  = "Adult patients with cancer, pooled across three previously conducted PK studies (Study 1 n=50, Study 2 n=7, Study 3 n=13).",
    dose_range     = "Oral sunitinib 25, 37.5, or 50 mg once daily; complete dosing histories were not available for each patient and pre-dose concentrations were handled via the Soy 'missing-dose' method during NONMEM estimation.",
    regions        = "Three multi-centre clinical-pharmacology studies (Netherlands Cancer Institute and partner centres).",
    notes          = "Pooled n = 1205 plasma samples (602 sunitinib + 603 SU12662) from 70 cancer patients. Baseline WT median 82 kg; per Yu 2015 Table 1 the WT range is 39-157 kg. Sex / race / age breakdown is not reported in the main text. Body weights for the simulation cohort follow a truncated log-normal distribution with mean 82.3 kg, SD 19.4 kg, truncated at 39-157 kg (Yu 2015 'Simulations of dosing regimens')."
  )

  ini({
    # ----------------------------------------------------------------------
    # Sunitinib (parent) PK -- Table 2, current-model column.
    # The model is 1-compartment for sunitinib with a hypothetical
    # pre-systemic / hepatic enzyme compartment in algebraic equilibrium
    # with the central compartment (Cliv = (ka*depot + Qh/Vc*central) /
    # (Qh + CL)). All sunitinib parameters are apparent (estimated
    # relative to the unknown oral bioavailability F).
    # ----------------------------------------------------------------------
    lka <- log(0.34);  label("Sunitinib absorption rate Ka (1/h)")                             # Table 2, Ka 0.34 (RSE 10.8%)
    lcl <- log(35.7);  label("Sunitinib apparent clearance CL/F at WT = 70 kg (L/h)")          # Table 2, CL 35.7 (RSE 5.7%)
    lvc <- log(1360);  label("Sunitinib apparent central volume Vc/F at WT = 70 kg (L)")       # Table 2, Vc 1360 (RSE 6.0%)

    # ----------------------------------------------------------------------
    # SU12662 (active metabolite) PK -- Table 2, current-model column.
    # 2-compartment disposition; apparent values (relative to the fixed
    # metabolite formation fraction fm = 0.21).
    # ----------------------------------------------------------------------
    lcl_su12662 <- log(17.1); label("SU12662 apparent clearance CL_M/(F * fm) at WT = 70 kg (L/h)")        # Table 2, CL 17.1 (RSE 7.4%)
    lvc_su12662 <- log(635);  label("SU12662 apparent central volume Vc_M/(F * fm) at WT = 70 kg (L)")     # Table 2, Vc 635 (RSE 13.1%)
    lq_su12662  <- log(20.1); label("SU12662 apparent inter-compartmental clearance Qi/(F * fm) at WT = 70 kg (L/h)") # Table 2, Qi 20.1 (RSE 32.6%)
    lvp_su12662 <- log(388);  label("SU12662 apparent peripheral volume Vp_M/(F * fm) at WT = 70 kg (L)")  # Table 2, Vp 388 (RSE 14.9%)

    # ----------------------------------------------------------------------
    # Allometric scaling -- Methods Eqs 3-4. Fixed exponents (no SE / CI
    # reported in the paper); reference WT = 70 kg.
    # ----------------------------------------------------------------------
    e_wt_cl <- fixed(0.75); label("Allometric exponent on CL, Qh, CL_SU12662, Qi_SU12662 (fixed; Methods Eq 3)")  # Methods Eq 3
    e_wt_vc <- fixed(1.0);  label("Allometric exponent on Vc, Vc_SU12662, Vp_SU12662 (fixed; Methods Eq 4)")       # Methods Eq 4

    # ----------------------------------------------------------------------
    # Inter-individual variability -- Table 2 + $OMEGA BLOCK(4) in the
    # NONMEM control stream (Appendix). The block carries non-zero
    # off-diagonals only between (Vc_sun, Vc_SU12662), (Vc_SU12662,
    # CL_SU12662), and (CL_SU12662, CL_sun); other pairs are estimated
    # to zero. Variances are computed from the published %CV via
    # omega^2 = log(CV^2 + 1); covariances from rho * sqrt(var_i * var_j)
    # using the correlations 0.48, 0.45, 0.53 in Table 2.
    # ETA order matches the NONMEM $OMEGA BLOCK(4):
    #   ETA1 = etalvc (Vc_sun),  ETA2 = etalvc_su12662 (Vc_SU12662),
    #   ETA3 = etalcl_su12662 (CL_SU12662), ETA4 = etalcl (CL_sun).
    # ----------------------------------------------------------------------
    etalvc + etalvc_su12662 + etalcl_su12662 + etalcl ~ c(
      0.0998,
      0.0816, 0.2891,
      0,      0.0978, 0.1632,
      0,      0,      0.0706, 0.1087
    ) # Table 2; variances via log(1 + CV^2) from %CV = 32.4 / 57.9 / 42.1 / 33.9; off-diagonals = rho * sqrt(var_i * var_j) with rho(Vc_sun,Vc_su)=0.48, rho(Vc_su,CL_su)=0.45, rho(CL_su,CL_sun)=0.53; other pairs fixed to 0 per OMEGA BLOCK(4) in the NONMEM Appendix

    # ----------------------------------------------------------------------
    # Residual error -- Table 2. Yu 2015 estimated study-specific sigma^2
    # values: sunitinib 0.06 (Study 1, n=50, 703 samples); 0.0188 init for
    # Studies 2+3 (NONMEM control stream); SU12662 0.03 (Study 1) and 0.01
    # (Studies 2+3). Per-study stratification does not generalise to a
    # downstream simulation, so a single proportional SD per output is
    # carried forward, using the Study 1 sigma^2 (the largest cohort).
    # propSd is the SD on the proportional term: SD = sqrt(sigma^2).
    # ----------------------------------------------------------------------
    propSd         <- sqrt(0.06); label("Proportional residual SD on sunitinib (Study 1; sqrt of Table 2 sigma^2 = 0.06)") # Table 2 (Study 1)
    propSd_su12662 <- sqrt(0.03); label("Proportional residual SD on SU12662  (Study 1; sqrt of Table 2 sigma^2 = 0.03)")  # Table 2 (Study 1)
  })

  model({
    # ------------------------------------------------------------------
    # Fixed structural constants from the paper (not estimated).
    # ------------------------------------------------------------------
    fm   <- 0.21    # Fraction of sunitinib converted to SU12662 (Methods + Table 2 "0.21 fix"; Houk 2009)
    Qh70 <- 80      # Hepatic blood flow at WT = 70 kg, L/h (Methods + Table 2 "80 fix")

    # ------------------------------------------------------------------
    # Allometric size factors (Methods Eqs 3-4; WT/70 baseline).
    # ------------------------------------------------------------------
    ascl <- (WT / 70)^e_wt_cl
    asv  <- (WT / 70)^e_wt_vc

    # ------------------------------------------------------------------
    # Individual PK parameters (typical * exp(eta) * allometric).
    # Sunitinib Ka has no IIV in the published OMEGA BLOCK.
    # ------------------------------------------------------------------
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * ascl
    vc <- exp(lvc + etalvc) * asv
    qh <- Qh70                  * ascl

    cl_su12662 <- exp(lcl_su12662 + etalcl_su12662) * ascl
    vc_su12662 <- exp(lvc_su12662 + etalvc_su12662) * asv
    q_su12662  <- exp(lq_su12662)                   * ascl
    vp_su12662 <- exp(lvp_su12662)                  * asv

    # ------------------------------------------------------------------
    # Micro-rate constants for the SU12662 sub-system (1/h).
    # ------------------------------------------------------------------
    kel_su <- cl_su12662 / vc_su12662
    k12_su <- q_su12662  / vc_su12662
    k21_su <- q_su12662  / vp_su12662

    # ------------------------------------------------------------------
    # Hypothetical hepatic enzyme-site concentration in quasi-steady-
    # state with the depot input and the sunitinib central compartment
    # via hepatic blood flow Qh (Yu 2015 Eq 6 + NONMEM $DES "Cliv =
    # (K12*A(1) + QH/V2*A(2)) / (QH + CLP)"). Sunitinib clearance acts
    # at the enzyme site; fraction fm of the cleared sunitinib is
    # converted to SU12662 and enters the SU12662 central compartment.
    # ------------------------------------------------------------------
    cliv <- (ka * depot + qh / vc * central) / (qh + cl)

    # ------------------------------------------------------------------
    # ODE system (Yu 2015 Eqs 7-10; NONMEM $DES in the Appendix).
    # Amounts are in mg, volumes in L, so plasma is in mg/L; multiplied
    # by 1000 to give ng/mL = ug/L.
    # ------------------------------------------------------------------
    d/dt(depot)               <- -ka * depot
    d/dt(central)             <-  qh * cliv - qh / vc * central
    d/dt(central_su12662)     <-  fm * cl * cliv -
                                  kel_su * central_su12662 -
                                  k12_su * central_su12662 +
                                  k21_su * peripheral1_su12662
    d/dt(peripheral1_su12662) <-  k12_su * central_su12662 -
                                  k21_su * peripheral1_su12662

    # ------------------------------------------------------------------
    # Observations (ng/mL = ug/L). The paper reports concentrations and
    # PK targets (50 ng/mL combined Cmin) in ng/mL throughout.
    # ------------------------------------------------------------------
    Cc         <- 1000 * central          / vc
    Cc_su12662 <- 1000 * central_su12662  / vc_su12662

    Cc         ~ prop(propSd)
    Cc_su12662 ~ prop(propSd_su12662)
  })
}
