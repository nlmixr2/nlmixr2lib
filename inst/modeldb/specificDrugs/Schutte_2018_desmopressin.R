Schutte_2018_desmopressin <- function() {
  description <- "Two-compartment apparent population PK model describing the time profile of endogenous factor VIII coagulant activity (FVIII:C) following a desmopressin (DDAVP) administration in nonsevere haemophilia A patients (Schutte 2018; final covariate model with FVIII-recent on baseline FVIII, V1 and CL). Desmopressin is the administered intervention; the apparent PK parameters describe the resulting endogenous FVIII:C release as if it were a unit-dose drug input (the source paper fixed the dose to unity because no FVIII concentrate was infused)."
  reference <- "Schutte LM, van Hest RM, Stoof SCM, Leebeek FWG, Cnossen MH, Kruip MJHA, Mathot RAA. Pharmacokinetic Modelling to Predict FVIII:C Response to Desmopressin and Its Reproducibility in Nonsevere Haemophilia A Patients. Thromb Haemost 2018;118(3):621-629. doi:10.1160/TH17-06-0390"
  vignette <- "Schutte_2018_desmopressin"
  units <- list(time = "h", dosing = "unit (DDAVP-triggered FVIII release; arbitrary)", concentration = "IU/mL")

  covariateData <- list(
    FVIIIRECENT = list(
      description        = "Most recently measured FVIII coagulant activity (FVIII:C) <=1 day before the desmopressin administration, in the absence of any treatment effect on the measurement. Per-occasion baseline-capacity index distinct from the model's observed FVIII:C time profile.",
      units              = "IU/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on baseline FVIII (+0.74), V1 (-0.61) and CL (-0.73) with reference 0.15 IU/mL (Schutte 2018 Table 2 and Eqs. 1-3; reference value is the study-population median FVIII-recent). Required input for simulation. The source paper also defines a binary missing-data branch with flat correction factors 1.2 / 1.1 / 0.78 on baseline FVIII / V1 / CL when FVIII-recent was unavailable for a fitted subject; this branch is documented in the validation vignette deviations and is NOT encoded in model() (see vignette).",
      source_name        = "FVIII-recent"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 128L,
    n_administrations = 142L,
    n_observations   = 623L,
    n_studies        = 1L,
    age_range        = "7-75 years",
    age_median       = "28 years",
    weight_range     = "26-120 kg",
    weight_median    = "75 kg",
    sex_female_pct   = 0,
    race_ethnicity   = "not reported in source (single-centre Dutch cohort)",
    disease_state    = "Nonsevere haemophilia A (FVIII activity >=0.01 IU/mL); mild HA 107/128 (84%) and moderate HA 21/128 (16%) with F8-gene mutations across A1 / A2 / A3 / B / C1 / C2 domains.",
    dose_range       = "0.3 ug/kg desmopressin intravenously over 30 minutes (126/137 administrations, 92%) or 300 ug intranasally (11/137 administrations, 8%).",
    regions          = "Netherlands (single-centre retrospective cohort at Erasmus University Medical Centre, Rotterdam).",
    indication       = "Desmopressin test dose 109/131 (83%), surgical prophylaxis 9/131 (7%), and treatment of bleeding 13/131 (10%).",
    notes            = "Haemophilia A is X-linked recessive; the registered cohort is essentially male, so sex_female_pct = 0. Baseline FVIII-recent (n = 120 with measurement available) had median 0.15 IU/mL (IQR 0.08 - 0.24); this is used as the reference value in the covariate equations. Sample times nominally 0.5, 1, 2, 3, 4, 6 and 24 h after the desmopressin administration (test-dose protocol); fewer time points in prophylaxis or bleed indications. 14 patients contributed two administrations on different occasions."
  )

  ini({
    # Structural parameters at the reference covariate value FVIIIRECENT = 0.15 IU/mL
    # (study-population median). Apparent (CL/F, V1/F, V2/F, Q/F) because the source
    # fixed the unit dose; the trigger drug is desmopressin and the modelled
    # response is endogenous FVIII:C.
    lka         <- log(3.8)    ; label("First-order absorption rate ka (1/h)")                                        # Schutte 2018 Table 2 final-model column: ka = 3.8, RSE 24%
    lcl         <- log(0.26)   ; label("Apparent clearance CL/F (L/h) at FVIIIRECENT = 0.15 IU/mL")                   # Schutte 2018 Table 2 final-model column: CL = 0.26 L/h, RSE 6.5%
    lvc         <- log(1.7)    ; label("Apparent central volume V1/F (L) at FVIIIRECENT = 0.15 IU/mL")                # Schutte 2018 Table 2 final-model column: V1 = 1.7 L, RSE 6.9%
    lq          <- log(0.11)   ; label("Apparent intercompartmental clearance Q/F (L/h)")                             # Schutte 2018 Table 2 final-model column: Q = 0.11 L/h, RSE 48%
    lvp         <- log(0.24)   ; label("Apparent peripheral volume V2/F (L)")                                         # Schutte 2018 Table 2 final-model column: V2 = 0.24 L, RSE 36%
    lbase_fviii <- log(0.15)   ; label("Typical baseline FVIII:C activity (IU/mL) at FVIIIRECENT = 0.15 IU/mL")       # Schutte 2018 Table 2 final-model column: Baseline FVIII = 0.15 IU/mL, RSE 4.6%

    # Covariate effects on baseline FVIII, V1 and CL (Schutte 2018 Eqs. 1-3 with
    # FVIIIRECENT available; the missing-FVIIIRECENT branch is documented in the
    # validation vignette and not encoded here).
    e_fviiirecent_base <-  0.74; label("Power exponent of FVIIIRECENT on baseline FVIII (unitless)")   # Schutte 2018 Table 2: FVIII recent on baseline = 0.74, RSE 18%
    e_fviiirecent_vc   <- -0.61; label("Power exponent of FVIIIRECENT on V1 (unitless)")               # Schutte 2018 Table 2: FVIII recent on V1      = -0.61, RSE 16%
    e_fviiirecent_cl   <- -0.73; label("Power exponent of FVIIIRECENT on CL (unitless)")               # Schutte 2018 Table 2: FVIII recent on CL      = -0.73, RSE 10%

    # Inter-individual variability. Schutte 2018 Table 2 final-model column reports
    # IIV as CV% (37% baseline FVIII, 43% V1, 50% CL) with a negative correlation
    # -0.87 between IIV(baseline FVIII) and IIV(V1). Variances are encoded via the
    # log-normal identity omega^2 = log(1 + CV^2):
    #   var_base = log(1 + 0.37^2) = log(1.1369) = 0.128306
    #   var_vc   = log(1 + 0.43^2) = log(1.1849) = 0.169666
    #   var_cl   = log(1 + 0.50^2) = log(1.25)   = 0.223144
    #   cov_base_vc = -0.87 * sqrt(0.128306 * 0.169666) = -0.128368
    # IOV (intraindividual variability) of 38% on baseline FVIII and 38% on V1
    # (also reported in Table 2, final-model column) is NOT implemented; see vignette
    # deviations.
    etalbase_fviii + etalvc ~ c(0.128306,
                                -0.128368, 0.169666)   # Schutte 2018 Table 2: IIV baseline 37% CV, IIV V1 43% CV, correlation -0.87
    etalcl                  ~ 0.223144                 # Schutte 2018 Table 2: IIV CL 50% CV

    # Residual error (final model with covariates, applied to the observed
    # FVIII:C = baseline + apparent increase).
    propSd <- 0.12 ; label("Proportional residual error (fraction of FVIII:C)")   # Schutte 2018 Table 2: proportional error = 12%, RSE 9.7%
    addSd  <- 0.018; label("Additive residual error (IU/mL)")                     # Schutte 2018 Table 2: additive error = 0.018 IU/mL, RSE 20%
  })
  model({
    # Individual PK parameters with the FVIIIRECENT power-scaling covariate
    # (Schutte 2018 Eqs. 1-3, FVIIIRECENT-available branch).
    base_fviii <- exp(lbase_fviii + etalbase_fviii) * (FVIIIRECENT / 0.15)^e_fviiirecent_base
    vc         <- exp(lvc         + etalvc)         * (FVIIIRECENT / 0.15)^e_fviiirecent_vc
    cl         <- exp(lcl         + etalcl)         * (FVIIIRECENT / 0.15)^e_fviiirecent_cl
    ka         <- exp(lka)
    vp         <- exp(lvp)
    q          <- exp(lq)

    # Micro-constants for the explicit two-compartment ODEs.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ODEs: desmopressin is administered as a unit perturbation into the depot;
    # the depot rate constant ka describes the time course of the desmopressin-
    # triggered FVIII:C release, which then distributes via the standard
    # two-compartment apparent disposition (Schutte 2018 Fig. 3).
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                k12 * central - k21 * peripheral1

    # Observation: total FVIII:C activity = baseline + DDAVP-triggered increase.
    # The dose is a unit perturbation, so central / vc carries the IU/mL scale
    # implied by the source's parameterisation (apparent V1/F in L).
    Cc <- base_fviii + central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
