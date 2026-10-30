LiricRajlic_2024_mRNA6231 <- function() {
  description <- "Preclinical (cynomolgus monkey). Semi-mechanistic kinetic-pharmacodynamic (KPD) model for mRNA-6231, a lipid-nanoparticle-encapsulated mRNA encoding a human interleukin-2 mutein fused to human serum albumin (HSA-IL2m), given subcutaneously. The lipid-nanoparticle and mRNA are below assay quantitation, so no drug concentration is fitted: a subcutaneous depot feeds a virtual transit compartment (unit volume) whose mRNA concentration drives saturable synthesis of HSA-IL2m protein into a central compartment. Protein in the central compartment equilibrates with an empirical effect compartment, and effect-site HSA-IL2m stimulates the production rate of a one-compartment indirect-response model for the percentage of regulatory T cells (FOXP3+ of CD4+). Structural parameters carry allometric body-weight scaling (reference 2.5 kg) so the monkey fit can be forward-projected to other body weights."
  reference <- paste(
    "Liric Rajlic I, Rong Y, Gaffney M, Faya P, Sun J, Liu C, Huang Y,",
    "Mahajan S, Hoffmann M, Smith M, Narasimhan B, Patel N.",
    "Translational kinetic-pharmacodynamics of mRNA-6231, an investigational",
    "mRNA therapeutic encoding mutein interleukin-2.",
    "CPT Pharmacometrics Syst Pharmacol. 2024;13(6):1067-1076.",
    "doi:10.1002/psp4.13142.",
    sep = " "
  )
  vignette <- "LiricRajlic_2024_mRNA6231"
  units <- list(
    time = "h",
    dosing = "ug",
    concentration = "ng/mL"
  )

  # `treg` holds the percentage of FOXP3+ regulatory T cells among CD4+ T
  # cells -- a dimensionless indirect-response state, not a mass compartment,
  # and with no canonical counterpart.
  paper_specific_compartments <- c("treg")

  compartmentData <- list(
    depot = list(
      analyte = "mRNA-6231",
      units = "ug",
      specimen = "administration site",
      verified = TRUE
    ),
    transit1 = list(
      analyte = "mRNA-6231",
      units = "ug",
      specimen = "not applicable",
      verified = TRUE
    ),
    central = list(
      analyte = "HSA-IL2m protein",
      units = "ng",
      specimen = "plasma",
      verified = TRUE
    ),
    effect = list(
      analyte = "HSA-IL2m protein (effect site)",
      units = "ng/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    treg = list(
      analyte = "regulatory T cells (% FOXP3+ of CD4+)",
      units = "%",
      specimen = "blood cell",
      verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric scaling of ka, CLt, Vt, V, CLIL2 and IL2kp (Table 2, equations 7-11), each normalised to the 2.5 kg reference body weight of the cynomolgus monkey analysis population. Fixed exponents: -0.25 (ka, IL2kp), 0.65 (CLt), 1.00 (Vt, V), 0.85 (CLIL2).",
      source_name = "WT"
    )
  )

  population <- list(
    species = "cynomolgus monkey (Macaca fascicularis)",
    n_subjects = 44L,
    n_studies = 3L,
    age_range = "2-5 years",
    weight_median = "2.5 kg (allometric reference)",
    sex_female_pct = 0,
    disease_state = "Healthy, naive and non-naive male cynomolgus monkeys",
    dose_range = "Single or repeat (Days 1, 15, 29) subcutaneous mRNA-6231. Study 1: 1, 3, 6 or 10 ug/kg on Days 1, 15 and 29 (N = 16). Study 2: single 3 or 10 ug/kg, or a flat 10 ug total dose (N = 12). Study 3: single 10 or 30 ug/kg, or a flat 100 ug total dose (N = 16).",
    notes = "Three non-GLP studies (Methods, Study design). HSA-IL2m protein was measured in plasma by ELISA; percentage Tregs (FOXP3+ of CD4+) was measured in cryopreserved PBMCs by flow cytometry. Systemic lipid-nanoparticle and mRNA concentrations were below the assay LLOQ and were not fitted, hence the kinetic-pharmacodynamic (dose-driven) rather than concentration-driven structure."
  )

  ini({
    # ---------------------------------------------------------------- PK ----
    # mRNA absorption and virtual-transit disposition (Table 2, equations 1-2).
    lka <- log(0.0749)
    label("Subcutaneous absorption rate ka (1/h)") # Table 3, ka = 0.0749 1/h (95% CI 0.06-0.08)
    lcl_transit <- log(12.865)
    label("Clearance of mRNA-6231 from the virtual transit compartment CLt (mL/h)") # Table 3, CL = 12.865 mL/h (95% CI 10.09-15.63)
    lvt <- fixed(log(1))
    label("Volume of the virtual transit compartment Vt (mL)") # Table 3, Vt = 1.0 mL (fixed)
    lvc <- log(171.56)
    label("Volume of distribution of HSA-IL2m protein V (mL)") # Table 3, V = 171.56 mL (95% CI 129.59-213.53)

    # Saturable synthesis of HSA-IL2m protein driven by transit-compartment
    # mRNA concentration, with first-order protein elimination (Table 2, eq. 3).
    lksyn <- log(4570.3)
    label("HSA-IL2m synthesis rate coefficient IL2kp (ng/h)") # Table 3, IL2kp = 4570.3 (95% CI 813.93-8326.6)
    lkm_syn <- log(0.1312)
    label("Transit-compartment mRNA concentration for half-maximal HSA-IL2m synthesis IL2k50 (ug/mL)") # Table 3, IL2k50 = 0.1312 ug/mL (95% CI 0.041-0.221)
    lcl <- log(1.529)
    label("Clearance of HSA-IL2m protein CLIL2 (mL/h)") # Table 3, CLIL2 = 1.529 mL/h (95% CI 1.23-1.91)

    # --------------------------------------------------------------- PD ----
    # Effect-compartment equilibration and indirect-response Treg turnover
    # (Table 2, equations 4-6).
    lke0 <- log(0.011)
    label("Effect-compartment equilibration rate ke0 (1/h)") # Table 3, ke0 = 0.011 1/h (95% CI 0.010-0.013)
    lkin <- log(0.0597)
    label("Treg production rate kin-Treg (%/h)") # Table 3, kin-Treg = 0.0597 %/h (95% CI 0.048-0.070)
    lkout <- log(0.0152)
    label("Treg loss rate kout-Treg (1/h)") # Table 3, kout-Treg = 0.0152 1/h (95% CI 0.013-0.017)
    lemax <- log(2.84)
    label("Maximal fold increase in Treg production Emax (unitless)") # Table 3, Emax = 2.84 (95% CI 2.38-4.31)
    lhill <- log(2.47)
    label("Sigmoidicity of the HSA-IL2m effect on Treg production gam (unitless)") # Table 3, gam = 2.47 (95% CI 2.04-3.15)
    lec50 <- log(5.08)
    label("Effect-site HSA-IL2m concentration for half-maximal Treg stimulation EC50 (ng/mL)") # Table 3, EC50 = 5.08 ng/mL (95% CI 4.61-6.12)

    # ---------------------------------------- Allometric scaling (fixed) ----
    # Fixed body-weight exponents, reference 2.5 kg (Table 2, equations 7-11;
    # text 'Extrapolation to humans': exponents fixed at 0.85, 0.65, 1 and
    # -0.25).
    e_wt_ka <- fixed(-0.25)
    label("Allometric exponent of body weight on ka (unitless)") # Table 2, equation 7
    e_wt_cl_transit <- fixed(0.65)
    label("Allometric exponent of body weight on CLt (unitless)") # Table 2, equation 9
    e_wt_vt <- fixed(1.00)
    label("Allometric exponent of body weight on Vt (unitless)") # Table 2, equation 10
    e_wt_vc <- fixed(1.00)
    label("Allometric exponent of body weight on V (unitless)") # Table 2, equation 10
    e_wt_cl <- fixed(0.85)
    label("Allometric exponent of body weight on CLIL2 (unitless)") # Table 2, equation 8
    e_wt_ksyn <- fixed(-0.25)
    label("Allometric exponent of body weight on IL2kp (unitless)") # Table 2, equation 11

    # ------------------------------------------------- Interindividual ------
    # Exponential IIV reported as %CV (Table 3). Variances use the exact
    # log-normal relation omega^2 = log(1 + CV^2).
    etalksyn ~ 0.233452 # Table 3, 'IIV on IL2kp' 51.3% (shrinkage 26.0) -> log(1 + 0.513^2)
    etalcl ~ 0.050676 # Table 3, 'IIV on CLIL2' 22.8% (shrinkage 31.0) -> log(1 + 0.228^2)
    etalkin ~ 0.125092 # Table 3, 'IIV on kin-Treg' 36.5% (shrinkage 4.0) -> log(1 + 0.365^2)
    etalemax ~ 0.087281 # Table 3, 'IIV on Emax' 30.2% (shrinkage 4.0) -> log(1 + 0.302^2)

    # ------------------------------------------------- Residual error -------
    propSd <- 0.434
    label("HSA-IL2m proportional residual error (fraction)") # Table 3, 'Proportional error IL-2' 43.4% CV
    propSd_treg <- 0.271
    label("Percentage-Treg proportional residual error (fraction)") # Table 3, 'Proportional error Treg' 27.1% CV
  })

  model({
    # 1. Allometric body-weight scaling, reference 2.5 kg (Table 2, eqs. 7-11)
    wtr <- WT / 2.5
    ka <- exp(lka) * wtr^e_wt_ka
    cl_transit <- exp(lcl_transit) * wtr^e_wt_cl_transit
    vt <- exp(lvt) * wtr^e_wt_vt
    vc <- exp(lvc) * wtr^e_wt_vc
    cl <- exp(lcl + etalcl) * wtr^e_wt_cl
    ksyn <- exp(lksyn + etalksyn) * wtr^e_wt_ksyn

    # 2. Weight-independent individual parameters
    km_syn <- exp(lkm_syn)
    ke0 <- exp(lke0)
    kin <- exp(lkin + etalkin)
    kout <- exp(lkout)
    emax <- exp(lemax + etalemax)
    hill <- exp(lhill)
    ec50 <- exp(lec50)

    # 3. mRNA disposition (Table 2, equations 1-2). Ct is the transit-
    # compartment mRNA concentration driving protein synthesis (equation 3).
    ct <- transit1 / vt
    d/dt(depot) <- -ka * depot
    d/dt(transit1) <- ka * depot - cl_transit * ct

    # 4. Saturable HSA-IL2m synthesis and first-order elimination (equation 3).
    # The CLIL2 allometric relation (equation 8) also carries a dose effect
    # (Dose/0.0045)^nonlin, but the 'nonlin' exponent is not reported anywhere
    # in the paper; the term is omitted (set to 1). With it at unity this model
    # reproduces the observed monkey NCA of Table 1. See vignette Errata.
    d/dt(central) <- ksyn * ct / (ct + km_syn) - cl * central / vc
    Cc <- central / vc

    # 5. Effect compartment equilibrating with plasma HSA-IL2m (equation 4).
    d/dt(effect) <- ke0 * (Cc - effect)

    # 6. Indirect-response Treg turnover, stimulated production (equation 5);
    # drug-free steady state kin/kout at baseline (equation 6).
    treg(0) <- kin / kout
    d/dt(treg) <- kin * (1 + emax * effect^hill / (ec50^hill + effect^hill)) -
      kout * treg

    # 7. Observations and residual error
    Cc ~ prop(propSd)
    treg ~ prop(propSd_treg)
  })
}
