Kloft_2004_sibrotuzumab <- function() {
  description <- "Two-compartment population PK model for sibrotuzumab in adults with metastatic FAP-positive cancer (Kloft 2004), with parallel linear and Michaelis-Menten elimination from the central compartment and a fixed linear body-weight covariate (centered at 75 kg) on linear CL, central and peripheral volumes, and Vmax."
  reference <- paste(
    "Kloft C, Graefe E-U, Tanswell P, Scott AM, Hofheinz R, Amelsberg A, Karlsson MO. (2004).",
    "Population pharmacokinetics of sibrotuzumab, a novel therapeutic monoclonal antibody,",
    "in cancer patients. Invest New Drugs 22(1):39-52.",
    "doi:10.1023/B:DRUG.0000006173.72210.1c.",
    "DDMORE Foundation Model Repository: DDMODEL00000195.",
    sep = " "
  )
  vignette <- "Kloft_2004_sibrotuzumab"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")
  ddmore_id    <- "DDMODEL00000195"
  replicate_of <- NULL

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear covariate centered at 75 kg, applied as `param * (1 + e_wt_<param> * (WT - 75))` on linear CL, central volume V1, peripheral volume V2, and Michaelis-Menten Vmax. Reference WT = 75 kg per the DDMORE Executable_sibrotuzumab.mdl GROUP_VARIABLES block. The four `BETA_*_WT` effect coefficients are declared `fix = true` in the .mdl PARAMETERS block (i.e., not estimated alongside the other thetas in this DDMORE encoding).",
      source_name        = "WT"
    )
  )

  population <- list(
    n_subjects     = NA_integer_,
    n_studies      = NA_integer_,
    age_range      = NA_character_,
    weight_range   = NA_character_,
    sex_female_pct = NA_real_,
    disease_state  = "Adults with metastatic fibroblast-activation-protein (FAP)-positive cancers (Kloft 2004 enrolled patients with metastatic colorectal, non-small-cell lung, and head-and-neck carcinomas in a Phase I/II program of repeated weekly IV sibrotuzumab infusions).",
    dose_range     = "Repeated weekly intravenous sibrotuzumab infusions (Phase I dose escalation per Kloft 2004). Detailed dose levels are not reproduced in the DDMORE bundle; the bundle's Simulated_sibrotuzumab.csv ships a single 80 mg / 1-hour IV infusion per subject as a smoke test.",
    regions        = NA_character_,
    notes          = "Population demographic detail (n_subjects, weight range, age, sex distribution) is not reproduced in the DDMORE Foundation Model Repository bundle for DDMODEL00000195, and the original Kloft 2004 publication is not on disk in this worktree. The bundle's Simulated_sibrotuzumab.csv carries 20 virtual subjects with WT in {70, 80, 90, 100} kg drawn evenly across IDs as a regression-test cohort. See the validation vignette's Errata section for the full list of bundle-versus-publication caveats."
  )

  ini({
    # Structural PK parameters - DDMODEL00000195 Executable_sibrotuzumab.mdl PARAMETERS / STRUCTURAL block
    # values, mirrored in the rendered NMTRAN $THETA in Output_simulated_nca_simulation.1.lst lines 72-83.
    # The bundle ships only a $SIMULATION listing (no $ESTIMATION run), so these are the publication-
    # derived point values that drive the bundle's simulation rather than re-fitted final estimates.
    lcl   <- log(0.0221) ; label("Linear clearance CL (L/h)")                        # POP_CLL  in .mdl PARAMETERS
    lvc   <- log(4.13)   ; label("Central volume of distribution V1 (L)")            # POP_V1   in .mdl PARAMETERS
    lvp   <- log(3.19)   ; label("Peripheral volume of distribution V2 (L)")         # POP_V2   in .mdl PARAMETERS
    lq    <- log(0.0376) ; label("Inter-compartmental clearance Q (L/h)")            # POP_Q    in .mdl PARAMETERS
    lvmax <- log(0.0338) ; label("Michaelis-Menten maximum elimination rate Vmax (mg/h)") # POP_Vmax in .mdl PARAMETERS
    lkm   <- log(8.0)    ; label("Michaelis-Menten constant Km (mg/L)")              # POP_Km   in .mdl PARAMETERS

    # Body-weight covariate effects (linear form, centered at 75 kg). All four are declared
    # `fix = true` in Executable_sibrotuzumab.mdl PARAMETERS, so they are fixed at the values
    # below rather than estimated.
    e_wt_cl   <- fixed(0.0182)   ; label("Linear WT effect on CL: cl   = cl   * (1 + e_wt_cl   * (WT - 75)) (1/kg)")  # BETA_CLL_WT  in .mdl PARAMETERS (fix = true)
    e_wt_vc   <- fixed(0.0125)   ; label("Linear WT effect on V1: vc   = vc   * (1 + e_wt_vc   * (WT - 75)) (1/kg)")  # BETA_V1_WT   in .mdl PARAMETERS (fix = true)
    e_wt_vp   <- fixed(0.0105)   ; label("Linear WT effect on V2: vp   = vp   * (1 + e_wt_vp   * (WT - 75)) (1/kg)")  # BETA_V2_WT   in .mdl PARAMETERS (fix = true)
    e_wt_vmax <- fixed(0.00934)  ; label("Linear WT effect on Vmax: vmax = vmax * (1 + e_wt_vmax * (WT - 75)) (1/kg)") # BETA_Vmax_WT in .mdl PARAMETERS (fix = true)

    # Inter-individual variability. The DDMORE bundle reports each PPV as a standard deviation
    # on the log-normal eta scale (`type is sd` in .mdl VARIABILITY); the corresponding
    # log-scale variance (omega^2) used by nlmixr2 is the square of the SD. The rendered
    # NMTRAN $OMEGA echoes this as `STANDARD` with the SD value (Output_simulated_*.lst lines
    # 84-87) and the simulation step's "INITIAL ESTIMATE OF OMEGA" reports the squared values
    # (0.3249, 0.04, 0.04, 0.0841).
    #   PPV_CLL   SD = 0.57  -> omega^2 = 0.3249
    #   PPV_V1    SD = 0.20  -> omega^2 = 0.04
    #   PPV_V2    SD = 0.20  -> omega^2 = 0.04
    #   PPV_Vmax  SD = 0.29  -> omega^2 = 0.0841   (FIXED in .mdl)
    etalcl   ~ 0.3249
    etalvc   ~ 0.04
    etalvp   ~ 0.04
    etalvmax ~ fixed(0.0841)

    # Combined residual error. The DDMORE encoding uses the MDL `combinedError1` form,
    # rendered to NMTRAN as `W = RUV_ADD + RUV_PROP*IPRED; Y = IPRED + W*EPS(1)` with
    # var(EPS) = 1 (Output_simulated_*.lst lines 67-68 and $SIGMA 1.0 FIX line 88).
    # That is the linear-SD combined-error parameterization (`combined1()` in nlmixr2),
    # not the default Pythagorean-SD combined2 form.
    addSd  <- 0.093  ; label("Additive residual error (mg/L)")          # RUV_ADD  in .mdl PARAMETERS
    propSd <- 0.0491 ; label("Proportional residual error (fraction)")  # RUV_PROP in .mdl PARAMETERS
  })
  model({
    # Individual PK parameters with the linear body-weight covariate centered at 75 kg.
    cl   <- exp(lcl   + etalcl)   * (1 + e_wt_cl   * (WT - 75))
    vc   <- exp(lvc   + etalvc)   * (1 + e_wt_vc   * (WT - 75))
    vp   <- exp(lvp   + etalvp)   * (1 + e_wt_vp   * (WT - 75))
    vmax <- exp(lvmax + etalvmax) * (1 + e_wt_vmax * (WT - 75))
    q    <- exp(lq)
    km   <- exp(lkm)

    # Concentrations in central / peripheral compartments. Dose units: mg; volumes: L;
    # so Cc and Cp are in mg/L (= ug/mL). IV dosing targets the central compartment
    # (the DDMORE bundle's NMTRAN data records use AMT + RATE on COMP1).
    Cc <- central / vc
    Cp <- peripheral1 / vp

    # Two-compartment IV PK with parallel linear (CL) and saturable Michaelis-Menten
    # (Vmax / Km) elimination from the central compartment. ODEs reproduce the
    # Executable_sibrotuzumab.mdl MODEL_PREDICTION DEQ block (and equivalently the
    # rendered $DES in Output_simulated_*.lst lines 53-59).
    d/dt(central)     <- q * (Cp - Cc) - cl * Cc - vmax * Cc / (km + Cc)
    d/dt(peripheral1) <- q * (Cc - Cp)

    Cc ~ add(addSd) + prop(propSd) + combined1()
  })
}
