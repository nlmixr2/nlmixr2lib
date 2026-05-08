Garmann_2017_BAY81_8973 <- function() {
  description <- "Two-compartment population PK model for BAY 81-8973 (Kovaltry, full-length unmodified recombinant human factor VIII) in patients with severe haemophilia A aged 1-61 years pooled from the LEOPOLD I, II and Kids trials (Garmann 2017). Final model uses NONMEM M3 likelihood for samples below the chromogenic-assay limit of quantitation (1.5 IU/dL)."
  reference <- "Garmann D, McLeay S, Shah A, Vis P, Maas Enriquez M, Ploeger BA. Population pharmacokinetic characterization of BAY 81-8973, a full-length recombinant factor VIII: lessons learned -- importance of including samples with factor VIII levels below the quantitation limit. Haemophilia. 2017 Jul;23(4):528-537. doi:10.1111/hae.13192. PMID:28419647."
  vignette <- "Garmann_2017_BAY81_8973"
  units <- list(time = "hour", dosing = "IU", concentration = "IU/dL")

  covariateData <- list(
    LBM = list(
      description        = "Lean body weight (canonical column LBM; source paper uses LBW)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on CL and Vc with reference LBW = 51.1 kg, the cohort median (Garmann 2017 Table 2 footnote; Table 1 reports median LBW 51.1 kg). Body-composition formula not stated in the paper; in haemophilia popPK literature LBW is most commonly computed via the Hume (1966) or James (1976) formula. Stored under canonical LBM (lean body mass) per inst/references/covariate-columns.md; LBW and LBM refer to the same quantity. Body weight, height, and BMI were screened during covariate analysis but only LBW was retained in the final model.",
      source_name        = "LBW"
    )
  )

  population <- list(
    n_subjects     = 183L,
    n_studies      = 3L,
    age_range      = "1-61 years",
    age_median     = "22 years (mean 23.3, SD 14.6)",
    weight_range   = "11-124 kg",
    weight_median  = "60 kg (mean 57.7, SD 24.8)",
    height_range   = "74-192 cm",
    height_median  = "170 cm (mean 160, SD 25.8)",
    bmi_range      = "13-38.3 kg/m^2",
    bmi_median     = "20.4 kg/m^2 (mean 21.3, SD 5.14)",
    lbw_range      = "9.25-79.2 kg",
    lbw_median     = "51.1 kg (mean 46.2, SD 16.7); used as the covariate-centering reference",
    sex_female_pct = 0,
    race_ethnicity = c(White = 132L, Asian = 31L, Black = 10L, Hispanic = 9L, Other = 1L),
    disease_state  = "Severe haemophilia A (FVIII activity < 1% by one-stage clotting assay) with no history of FVIII inhibitors and >= 150 (LEOPOLD I and II) or >= 50 (LEOPOLD Kids part A) prior FVIII exposure days",
    dose_range     = "Single and repeat IV infusions over ~10 minutes; simulations in the paper used 25 or 50 IU/kg twice weekly",
    regions        = "Multinational LEOPOLD trials (LEOPOLD I [adults/adolescents 12-61 y], LEOPOLD II [adults/adolescents 12-61 y], LEOPOLD Kids part A [<= 12 y])",
    notes          = "Pooled analysis of 1535 chromogenic FVIII activity observations from 183 male haemophilia A patients across the 3 LEOPOLD trials; 16.5% of samples were below the lower limit of quantitation (1.5 IU/dL for the majority; 3 IU/dL for a small number). Sex is essentially all-male because haemophilia A is X-linked. Haemophilia B is a different disease (factor IX deficiency); see Koopman_2023_factorix.R for the analogous FIX model."
  )

  ini({
    # Structural parameters - typical values for the paper's reference patient
    # (LBW = 51.1 kg, the cohort median). Units: CL and CLp in dL/h, Vc and Vp
    # in dL. Source: Garmann 2017 Table 2 ("final BAY 81-8973 model" column).
    # Garmann 2017 Methods reports a log-normal BSV model
    # (P_j = P_pop * exp(eta_j)), so structural fixed effects are stored on the
    # log scale (lcl / lvc / lq / lvp).
    lcl <- log(1.88); label("Clearance for the reference 51.1 kg LBW patient (CL, dL/h)") # Garmann 2017 Table 2: CL = 1.88 dL/h
    lvc <- log(30.0); label("Central volume of distribution for the reference 51.1 kg LBW patient (Vc, dL)") # Garmann 2017 Table 2: Vc = 30.0 dL
    lq  <- log(1.90); label("Intercompartmental clearance (CLp, dL/h)") # Garmann 2017 Table 2: CLp = 1.90 dL/h
    lvp <- log(6.37); label("Peripheral volume of distribution (Vp, dL)") # Garmann 2017 Table 2: Vp = 6.37 dL

    # Covariate effects: power-form (allometric-style) effect of LBW on CL and Vc,
    # both centered at the median LBW of 51.1 kg. Form (Garmann 2017 Table 2
    # footnote): TV_par = theta_par * (LBW / 51.1)^theta_xx.
    e_lbw_cl <- 0.610; label("Power exponent of LBW on CL (unitless)") # Garmann 2017 Table 2: 0.610 (95% CI 0.45-0.75)
    e_lbw_vc <- 0.950; label("Power exponent of LBW on Vc (unitless)") # Garmann 2017 Table 2: 0.950 (95% CI 0.89-1.02)

    # Inter-individual variability. Garmann 2017 reports BSV as %CV with the
    # log-normal model P_j = P_pop * exp(eta_j); for log-normal, the variance on
    # the eta scale is omega^2 = log((CV/100)^2 + 1).
    #   CV(CL) = 37.0% -> omega^2 = log(1 + 0.370^2) = 0.12826
    #   CV(Vc) = 11.2% -> omega^2 = log(1 + 0.112^2) = 0.01247
    # Garmann 2017 tested a CL:Vc covariance during model development
    # ("Covariate analysis and final model") but did not retain it
    # (deltaOBJ = 4.8, below the 7.9 retention threshold), so etalcl and etalvc
    # are independent. No IIV was reported on CLp (Q) or Vp.
    etalcl ~ 0.12826 # Garmann 2017 Table 2: BSV CL = 37.0% CV
    etalvc ~ 0.01247 # Garmann 2017 Table 2: BSV Vc = 11.2% CV

    # Combined proportional + additive residual unexplained variability
    # (Garmann 2017 Methods, Eq. for C_ij): C_ij = C_hat * (1 + eps_prop) + eps_add.
    # propSd is the SD of eps_prop (CV expressed as fraction); addSd is the SD of
    # eps_add in IU/dL.
    propSd <- 0.267; label("Proportional residual error (fraction)") # Garmann 2017 Table 2: proportional RUV = 26.7% CV
    addSd  <- 1.10;  label("Additive residual error (IU/dL)") # Garmann 2017 Table 2: additive RUV = 1.10 IU/dL
  })
  model({
    # LBW power scaling on CL and Vc with reference 51.1 kg.
    cl <- exp(lcl + etalcl) * (LBM / 51.1)^e_lbw_cl
    vc <- exp(lvc + etalvc) * (LBM / 51.1)^e_lbw_vc
    q  <- exp(lq)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment IV model: dose enters central directly. BAY 81-8973 is
    # given as a short IV infusion (the paper fixed duration to 0.167 h = 10 min
    # when missing from dosing records); the infusion duration is supplied
    # through the user's event table, not the structural model.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # FVIII activity: dose in IU and central volume in dL -> central / vc has
    # units IU/dL (chromogenic FVIII activity, equivalent to % of normal FVIII).
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
