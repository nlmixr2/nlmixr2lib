Berthaud_2025_cefazolin <- function() {
  description <- "One-compartment population PK model for intravenous cefazolin in children undergoing maintenance hemodialysis for kidney failure (Berthaud 2025, n = 6 patients aged 1.3-14.6 years weighing 11.4-51 kg, 83 total plasma concentrations). Total apparent clearance is the SUM of a residual (non-renal plus residual renal) elimination clearance and a dialysis clearance that is switched on only while a hemodialysis session is running, gated by the time-varying RRT_HEMODIAL_ACTIVE covariate. Body weight enters both the residual clearance and the volume of distribution allometrically, normalized to 70 kg, with the exponents FIXED at 0.75 and 1. Dialysis membrane surface area (FILT_SA) drives the dialysis-clearance arm through an estimated power function centred at 1 m2 and explained essentially all of the between-subject variability on that arm, so no IIV is carried there; population dialysis clearance is more than 10-fold the population residual clearance. Serum albumin, fat-free mass, age, blood flow rate, ultrafiltration volume, vascular access type and RRT technique were tested and not retained. Estimated in Monolix 2023R1 by SAEM."
  reference <- "Berthaud R, Urien S, Krid S, Foissac F, Oualha M, Thy M, Boyer O, Beranger A, Hirt D, Benaboud S, Treluyer JM, Bouazza N. Cefazolin population pharmacokinetics in children undergoing maintenance hemodialysis for kidney failure. Antimicrob Agents Chemother. 2025;69(11):e00451-25. doi:10.1128/aac.00451-25. PMCID PMC12587577. ClinicalTrials.gov NCT02539407 (Optimome study). All parameter estimates from Table 2; the structural equations from the Results 'Population PK modeling' display-equation block."
  vignette <- "Berthaud_2025_cefazolin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Berthaud 2025 Methods ("Cefazolin and albumin assays"):
  # cefazolin TOTAL plasma concentrations were quantified by HPLC-UV,
  # calibration range 0.5-200 mg/L. The model carries no protein-binding
  # sub-model -- serum albumin was tested as a covariate and not retained --
  # so the single state holds total cefazolin.
  compartmentData <- list(
    central = list(analyte = "cefazolin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Berthaud 2025 Table 1: 11.4-51 kg across the six patients (11.4, 12.3, 23.2, 24.5, 47.5, 51.0). Enters both the residual elimination clearance and the volume of distribution as an allometric power ratio normalized to 70 kg, with the exponents FIXED (not estimated) at 0.75 for CL and 1 for Vd (Methods: 'The effect of BW was assessed according to the allometric rule with beta values fixed at 0.75 for CL and CLdial and 1 for Vd'; Table 2 reports 'Fixed' in the RSE column for both). Note the 70 kg reference is an allometric standardization convention, NOT a cohort central value -- every patient in this paediatric cohort is far below it, so the reported CLpop and Vdpop are extrapolated typical values for a hypothetical 70 kg subject rather than values any subject exhibits. Adding the allometric BW scaling decreased the BIC by 14 units and reduced the BSV on CL, CLdial and Vd from 1.42 to 1.08, 0.95 to 0.62, and 0.51 to 0.2 respectively. Weight-based allometric scaling was ALSO tested on the dialysis clearance arm but became non-significant once dialysis membrane surface area was added, and was not retained in the final model. Fat-free mass was tested in place of total body weight and did not improve the fit.",
      source_name        = "BW"
    ),
    FILT_SA = list(
      description        = "Dialysis membrane (dialyser) surface area of the filter in use",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Berthaud 2025 Table 1 records a per-patient dialysis membrane surface area (DMSA) of 0.2, 0.3, 1, 1, 1.5 and 1.4/1.7 m^2 -- six distinct values used as a genuinely CONTINUOUS covariate in a power model, which is why this model founds the continuous FILT_SA canonical rather than reusing the discrete FILT_SA_MED / FILT_SA_LARGE indicators (the FILT_SA_MED register entry anticipates exactly this case). Membranes used were Sureflux 30L, Elisio 15H and Elisio 17H (Nipro) and FXpaed, Fx50 and Fx60 (Fresenius). Enters the dialysis-clearance arm as an estimated power ratio centred at 1 m^2. Adding DMSA to CLdial made allometric BW on CLdial non-significant, decreased the BIC by a further 6 units, and reduced the BSV on CLdial towards 0 -- the paper concludes DMSA 'explained most of the interindividual variability on that parameter', which is why the final model carries NO IIV on the dialysis arm. Patient 6 switched membranes mid-study (1.4 then 1.7 m^2), so the covariate is time-varying within subject in the source data set. Meaningful only while RRT_HEMODIAL_ACTIVE = 1; the whole dialysis arm is gated off otherwise. For the paper's OWN dosing-regimen simulations DMSA was not measured but approximated from body weight as DMSA = 0.85 * BSA with BSA = (4 * BW + 7) / (BW + 90) (Results, 'Dosing regimen simulations'); the Discussion notes 0.85 sits mid-range of the recommended 75-100% of BSA and that CLdial would change proportionally under a different coefficient. That approximation is a simulation convenience, NOT the fitted data -- the Table 1 DMSA values do not equal 0.85 * BSA (e.g. patient 1, 23.2 kg: 0.85 * BSA = 0.75 vs an actual 1 m^2) -- so it is reproduced in the validation vignette's cohort construction, not in model().",
      source_name        = "DMSA"
    ),
    RRT_HEMODIAL_ACTIVE = list(
      description        = "Hemodialysis-active indicator (1 while a session is running, 0 in the interdialytic interval)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (interdialytic / no dialysis running)",
      notes              = "Time-varying within subject. Gates the ADDITIVE dialysis-clearance arm cl_hemodialysis so that total clearance reduces to the residual elimination clearance between sessions. The additive (rather than replacement) composition is settled by the paper's own printed final equation, 'CLtot_i = CL_i + CLdial_i' (Results, 'Population PK modeling'), and by the Introduction: 'during hemodialysis sessions, a dialysis clearance, itself variable, is ADDED to the patient residual elimination clearance (non-renal + residual renal clearance) of cefazolin'. This matters because the two sibling cefazolin/ceftriaxone hemodialysis models in the library (Duke_2024_cefazolin.R, Tsai_2023_ceftriaxone.R) use the REPLACEMENT rule instead, where the dialysis estimate is the total clearance during a session; reading Berthaud's additive arm as a replacement (or vice versa) would misstate clearance by the body baseline. The on/off indicator does NOT appear in the printed equation block, which states only the two clearance components and their sum; it is established by the surrounding text -- the Introduction sentence above, the Discussion's contrast between an 'interdialytic cefazolin half-life' and a 'dialysis half-life' (two distinct half-lives require two distinct clearance states), and the Figure 4 dosing simulations built on 'three dialysis sessions per week (days 3, 5 and 8)'. Session timing per patient is not tabulated in the paper; the cohort had a median 8 [4-12] (1-13) dialysis sessions each. Dialysis machines were Gambro AK200 Ultra S, Nikkiso DBB-05 and Fresenius 5008S CorDiax; techniques were hemodialysis (all patients) plus post-dilution hemofiltration in patient 6. Blood flow rate (median 142 mL/min) and ultrafiltration volume (median 550 mL) were tested as covariates on the dialysis arm and not retained.",
      source_name        = "(not a named source column; established from the study design)"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 6L,
    n_studies        = 1L,
    age_range        = "1.3-14.6 years (individual ages 1.3, 4.5, 11.0, 11.2, 14.3, 14.6; Table 1)",
    weight_range     = "11.4-51 kg (individual weights 11.4, 12.3, 23.2, 24.5, 47.5, 51.0; Table 1)",
    sex_female_pct   = 16.7,
    race_ethnicity   = "Not reported (single-centre paediatric nephrology department, Necker-Enfants Malades university hospital, Paris)",
    disease_state    = "Children with kidney failure established on maintenance hemodialysis for more than 6 months, all with methicillin-susceptible Staphylococcus aureus (MSSA) bloodstream infection. Causes of kidney failure: hemolytic uremic syndrome (n = 2), Denys-Drash syndrome (n = 1), Schimke immuno-osseous dysplasia (n = 1), unknown (n = 2). Patient 4 had undergone bilateral nephrectomy and was anephric, with an individual residual clearance of 0.0073 L/h -- almost five times lower than the second-lowest in the cohort -- and is visibly the least well described subject in the Figure 1 individual fits. All patients recovered with sterilized blood cultures.",
    dose_range       = "Cefazolin (1 g powder, Mylan) reconstituted to 20 mg/mL and given as 30-60 min intravenous intermittent infusions every 6, 8, 12, 24 or 48 h by programmable syringe pump. Median dose per infusion ranged 4.9-26.3 mg/kg across patients (per-patient means 8, 11.4, 12.6, 15.4, 19.7, 25.9 mg/kg; Table 1). Regimens were set by local experience and adjusted during treatment on the basis of plasma assay results.",
    regions          = "France (single centre, Paris)",
    renal_function   = "Kidney failure with GFR below 15 mL/min/1.73 m^2 (0.9 L/h/1.73 m^2) in all patients; GFR could not be precisely quantified at the time of the study. Residual (non-dialysis) elimination clearance nonetheless remained a material determinant of exposure -- a central finding of the paper.",
    notes            = "Prospective study conducted January 2018 - December 2019 as part of the wider Optimome study (ClinicalTrials.gov NCT02539407). Baseline characteristics from Table 1. 85 samples were available from 6 patients; 2 were discarded as outliers, leaving 83 analysed, with no concentrations below the limit of quantification. Bioanalysis: validated HPLC with UV detection at 280 nm (Kinetex C18 column, doripenem internal standard), calibration range 0.5-200 mg/L, inter- and intra-assay precision and accuracy < 15%. Maintenance hemodialysis modalities per Table 1: median 8 [4-12] dialysis sessions per patient, median blood flow rate 142 [110-260] mL/min, median ultrafiltration volume 550 [375-1200] mL. Model evaluation used prediction-corrected VPC and NPDE metrics from 500 Monte Carlo simulations per patient."
  )

  ini({
    # Structural parameters -- Berthaud 2025 Table 2, standardized to a body
    # weight of 70 kg. The published final equations (Results, 'Population PK
    # modeling') are:
    #   Vd_i(L)     = 14.6  * (BW_i / 70)
    #   CL_i(L/h)   = 0.186 * (BW_i / 70)^0.75
    #   CLdial_i(L/h) = 1.98 * (DMSA_i / 1)^1.26
    #   CLtot_i     = CL_i + CLdial_i
    lcl              <- log(0.186); label("Residual (non-renal + residual renal) elimination clearance at 70 kg (L/h)")  # Berthaud 2025 Table 2: CLpop = 0.186 L/h (RSE 45.7%)
    lcl_hemodialysis <- log(1.98);  label("Dialysis clearance at a 1 m^2 membrane surface area (L/h)")                    # Berthaud 2025 Table 2: CLdialpop = 1.98 L/h (RSE 3.74%)
    lvc              <- log(14.6);  label("Volume of distribution at 70 kg (L)")                                          # Berthaud 2025 Table 2: Vdpop = 14.6 L (RSE 10.9%)

    # Allometric exponents on body weight. Both are FIXED, not estimated:
    # Table 2 prints "Fixed" in the RSE column for each, and the Methods state
    # "The effect of BW was assessed according to the allometric rule with beta
    # values fixed at 0.75 for CL and CLdial and 1 for Vd".
    e_wt_cl <- fixed(0.75); label("Allometric exponent on (WT/70) for residual elimination clearance (unitless)")  # Berthaud 2025 Table 2: beta BW/CL = 0.75, Fixed
    e_wt_vc <- fixed(1);    label("Allometric exponent on (WT/70) for volume of distribution (unitless)")          # Berthaud 2025 Table 2: beta BW/Vd = 1, Fixed

    # Dialysis membrane surface area on the dialysis-clearance arm. Estimated
    # (RSE 1.06%), unlike the two allometric exponents above. Allometric BW on
    # CLdial was tested and became non-significant once this term was added.
    e_filt_sa_cl_hemodialysis <- 1.26; label("Power exponent on (FILT_SA/1) for dialysis clearance (unitless)")  # Berthaud 2025 Table 2: beta DMSA/CLdial = 1.26 (RSE 1.06%)

    # Between-subject variability. Table 2's footnote defines the reported
    # quantity explicitly: "omega, between-subject variability - square root of
    # the between-subjects variance omega^2", so the printed values are SDs on
    # the log scale and the variances below are their squares. The random
    # effects were "ascribed to an exponential distribution" (Methods).
    #   CL : omega = 1.07  -> omega^2 = 1.1449
    #   Vd : omega = 0.197 -> omega^2 = 0.038809
    # No IIV is carried on the dialysis arm: adding DMSA "reduced the BSV on
    # CLdial towards 0" and the final model reports no omega for it (Table 2).
    etalcl ~ 1.1449    # Berthaud 2025 Table 2: omegaCL = 1.07 (RSE 33.7%, shrinkage 8.43%)
    etalvc ~ 0.038809  # Berthaud 2025 Table 2: omegaVd = 0.197 (RSE 49.9%, shrinkage 20.6%)

    # Residual error: proportional only (Results: "A proportional model was
    # used to describe the residual variability"). Monolix's proportional error
    # is y = f * (1 + b * eps) with eps ~ N(0, 1), so the reported b maps
    # directly onto nlmixr2's prop() SD without transformation.
    propSd <- 0.398; label("Proportional residual error (fraction)")  # Berthaud 2025 Table 2: b (proportional) = 0.398 (RSE 9.33%)
  })

  model({
    # 1. Volume of distribution. Body weight enters allometrically, centred at
    #    the 70 kg allometric standard; the /70 and /1 denominators throughout
    #    are the paper's own printed centring values and are kept verbatim.
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc  # L

    # 2. Dialysis-clearance arm. Driven entirely by membrane surface area and
    #    carrying no random effect, because DMSA absorbed essentially all of
    #    the between-subject variability on this parameter.
    cl_hemodialysis <- exp(lcl_hemodialysis) * (FILT_SA / 1)^e_filt_sa_cl_hemodialysis  # L/h

    # 3. TOTAL apparent clearance. ADDITIVE composition, per the paper's own
    #    final equation CLtot_i = CL_i + CLdial_i, with the dialysis arm
    #    contributing only while a session is running. Interdialytically
    #    RRT_HEMODIAL_ACTIVE = 0 and cl collapses to the residual elimination
    #    clearance alone, which is what gives the model its two distinct
    #    half-lives (the Discussion contrasts an "interdialytic cefazolin
    #    half-life" with a "dialysis half-life").
    #
    #    IMPORTANT -- the gated sum MUST be assigned to `cl` itself, not to a
    #    separate `cl_total`. rxode2 recognises the joint presence of `cl` and
    #    `vc` and solves the one-compartment system analytically from that
    #    pair, discarding the explicit d/dt() right-hand side. Writing
    #    `cl <- <residual arm>` and then eliminating with `cl_total / vc`
    #    yields a model whose dialysis arm is silently INERT: the reported
    #    `cl_total` and `kel` columns look correct while the simulated
    #    concentrations decay at the interdialytic rate in both states.
    #    Verified in the validation vignette by a two-state gate check.
    #    Same hazard and same fix as Lee_2024_gentamicin_teigen.R.
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl +
      RRT_HEMODIAL_ACTIVE * cl_hemodialysis  # L/h

    kel <- cl / vc

    # 4. One-compartment disposition with intravenous input; no depot. Dose is
    #    delivered into 'central' with the duration set by the dosing event
    #    (30-60 min infusions in the source study).
    d/dt(central) <- -kel * central

    # 5. Observation. Dose in mg and vc in L give central/vc in mg/L, matching
    #    the total plasma cefazolin concentrations reported throughout.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
