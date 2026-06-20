Molto_2016_atazanavir_ritonavir <- function() {
  description <- paste(
    "Simultaneous one-compartment popPK model for oral atazanavir (ATV,",
    "parent / substrate) and ritonavir (RTV, sibling-drug suffix _rtv) in",
    "83 HIV-1-infected Caucasian adults receiving either ATV 400 mg or",
    "ATV 300 mg + RTV 100 mg once daily. Both drugs use a Savic transit-",
    "compartment absorption chain (ATV: N = 7, MTT = 0.80 h, ka = 2.05",
    "1/h; RTV: N = 11, MTT = 0.522 h, ka = 1.21 1/h) feeding a depot,",
    "followed by first-order elimination from a one-compartment central.",
    "ATV apparent clearance is exponentially inhibited by RTV plasma",
    "concentration: CL/F_ATV(t) = exp(lcl) * exp(-e_crtv_cl * C_RTV(t))",
    "with the unboosted CL/F_ATV = 11.7 L/h and inhibition coefficient",
    "0.296 L/mg. This functional form reproduces the paper's reported",
    "~18% reduction in ATV CL at the cohort-mean RTV concentration of",
    "0.63 mg/L and explains 17.5% of inter-individual variability in ATV",
    "CL. Demographic covariates (weight allometric, gender, age, TDF,",
    "HCV, dose-timing, AAG, albumin) were screened by GAM and tested in",
    "NONMEM but not retained; an Emax-form and a linear-form inhibition",
    "were also tested and rejected (Emax: unrealistic estimates; linear:",
    "biased fit). IIV on ka / CL/F / V/F is reported for both drugs with",
    "unusually large IIV on absorption (~200% CV) confirmed in the paper",
    "Results. ATV residual error is combined (27.0% proportional + 0.07",
    "mg/L additive); RTV residual error is proportional only (28.0%; the",
    "additive component of the initial combined error was deleted as",
    "negligible) (Molto 2016)."
  )
  reference <- paste(
    "Molto J, Estevez JA, Miranda C, Cedeno S, Clotet B, Valle M.",
    "Population pharmacokinetic modelling of the changes in atazanavir",
    "plasma clearance caused by ritonavir plasma concentrations in",
    "HIV-1 infected patients. Br J Clin Pharmacol.",
    "doi:10.1111/bcp.13072."
  )
  vignette <- "Molto_2016_atazanavir_ritonavir"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list()

  population <- list(
    species          = "human",
    n_subjects       = 83L,
    n_studies        = 1L,
    age_range        = "26-75 years",
    age_median       = "42 years",
    weight_range     = "40-91 kg",
    weight_median    = "70 kg",
    sex_female_pct   = 31.33,
    race_ethnicity   = "Caucasian (Spanish single-centre cohort)",
    disease_state    = paste(
      "HIV-1 infection on stable antiretroviral therapy for >=4 weeks",
      "(steady-state conditions). 30/83 with HCV coinfection; 5/83 with",
      "advanced liver fibrosis. Median AST 27 U/L (range 11-185);",
      "median ALT 30 U/L (range 7-279); median plasma proteins 7.4 g/dL",
      "(range 6.3-8.7); median albumin 4.3 g/dL (range 2.8-5.25);",
      "median alpha-1-acid glycoprotein 79 mg/dL (range 36-156)."
    ),
    dose_range       = paste(
      "ATV 400 mg orally once daily (unboosted, 30/83 patients) or",
      "ATV 300 mg + RTV 100 mg orally once daily (boosted, 53/83",
      "patients). Six dosing scenarios were simulated by the paper",
      "(ATV 400 QD; ATV 300/RTV 100 QD; ATV 300/RTV 50 QD; ATV 200/RTV",
      "100 QD; ATV 300 BID; ATV 200 BID)."
    ),
    regions          = "Spain (Catalonia, Hospital Universitari Germans Trias i Pujol, Badalona)",
    notes            = paste(
      "Cross-sectional study enrolled May 2004 to May 2009. 38/83",
      "patients were also taking tenofovir 300 mg once daily (2",
      "unboosted, 36 boosted). 16/83 took the medication at night",
      "(unwitnessed dose); 67/83 took it in the morning (witnessed",
      "dose). Per Molto 2016 Methods, neither dose-timing nor any",
      "other tested demographic / clinical covariate was retained in",
      "the final model. ATV LLOQ 0.044 mg/L, RTV LLOQ 0.05 mg/L",
      "(HPLC-PDA). <1% of ATV and ~5% of RTV samples were below LLOQ",
      "and excluded per Bergstrand 2009 (Methods ref 25). Median",
      "self-reported adherence 98.95%. Sampling: a baseline pre-dose",
      "sample plus 0.5, 1, 2, 4, 6, 8, 10, 12 h post-dose for",
      "morning-dosed patients (sparser, 12-22 h post-dose, for",
      "night-dosed patients)."
    )
  )

  ini({
    # ===================================================================
    # ATAZANAVIR (parent / substrate) structural parameters
    # ----- Molto 2016 Table 2, "Final model" columns for atazanavir -----
    # Both ATV and RTV use a Savic transit-compartment absorption chain
    # (Methods: "in order to describe the absorption profile, models ...
    #  with transit compartments as a more flexible function to describe
    #  the delay in the absorption were used").
    # ===================================================================
    lka  <- log(2.05)
    label("Atazanavir absorption rate constant ka (1/h)")             # Table 2 final ka_ATV = 2.05 1/h (RSE 19.5%)
    lmtt <- log(0.80)
    label("Atazanavir mean transit time MTT (h)")                      # Table 2 final MTT_ATV = 0.80 h (RSE 6.2%)
    lnn  <- fixed(log(7))
    label("Atazanavir number of transit compartments N (unitless, FIXED to integer 7)") # Table 2 N_transit_ATV = 7; no RSE reported because N was held at an integer selected by model search
    lcl  <- log(11.7)
    label("Atazanavir baseline apparent clearance in the absence of ritonavir, CL/F_ATV (L/h)") # Table 2 final CL/F_ATV = 11.7 L/h (RSE 6.8%); represents the unboosted typical value (interaction factor evaluates to 1 at C_RTV = 0)
    lvc  <- log(95.7)
    label("Atazanavir apparent volume of distribution, V/F_ATV (L)")  # Table 2 final V/F_ATV = 95.7 L (RSE 6.5%)

    # ===================================================================
    # RITONAVIR (sibling-drug suffix _rtv) structural parameters
    # ----- Molto 2016 Table 2, "Final model" columns for ritonavir -----
    # ===================================================================
    lka_rtv  <- log(1.21)
    label("Ritonavir absorption rate constant ka (1/h)")                # Table 2 final ka_RTV = 1.21 1/h (RSE 24.7%)
    lmtt_rtv <- log(0.522)
    label("Ritonavir mean transit time MTT (h)")                        # Table 2 final MTT_RTV = 0.522 h (RSE 3.8%)
    lnn_rtv  <- fixed(log(11))
    label("Ritonavir number of transit compartments N (unitless, FIXED to integer 11)") # Table 2 N_transit_RTV = 11; no RSE reported because N was held at an integer selected by model search
    lcl_rtv  <- log(9.68)
    label("Ritonavir apparent clearance, CL/F_RTV (L/h)")               # Table 2 final CL/F_RTV = 9.68 L/h (RSE 3.0%)
    lvc_rtv  <- log(70.5)
    label("Ritonavir apparent volume of distribution, V/F_RTV (L)")     # Table 2 final V/F_RTV = 70.5 L (RSE 8.8%)

    # ===================================================================
    # Drug-drug interaction: exponential inhibition of ATV apparent
    # clearance by RTV plasma concentration (Molto 2016 Methods,
    # "Final model" paragraph; Results: "the one incorporating the
    #  inhibition of CL/F ATV by RTV according to an exponential model
    #  ... best described the data"):
    #   CL/F_ATV(t) = exp(lcl) * exp(-e_crtv_cl * C_RTV(t))
    # Reproduces ~17-18% reduction in ATV CL/F at the cohort-mean RTV
    # plasma concentration of 0.63 mg/L (Molto 2016 Results: "A mean
    # RTV plasma concentration of 0.63 mg/L predicted an 18% decrease
    # in ATV clearance"). Linear and Imax inhibition forms were tested
    # and rejected (Methods + Results). At C_RTV = 0 the factor
    # evaluates to 1 so the unboosted CL/F_ATV equals exp(lcl) = 11.7
    # L/h directly.
    # ===================================================================
    e_crtv_cl <- 0.296
    label("Exponential inhibition coefficient of RTV plasma concentration on ATV CL/F (L/mg)") # Table 2 final theta_C(RTV),CL/F = 0.296 (RSE 3.5%)

    # ===================================================================
    # IIV - log-normal; omega^2 = log(1 + CV^2). Reported CV% from
    # Molto 2016 Table 2 "Final model" columns. Note the unusually
    # large IIV on absorption (~200% CV for both drugs) is reported
    # in the source and confirmed in Results ("The highest
    # interindividual variability estimated in our model was
    # associated with the ATV ka").
    # ===================================================================
    etalka     ~ log(1 + 2.002^2)  # IIV ka_ATV = 200.2% (RSE 22.4%);   omega^2 = log(1 + 2.002^2)
    etalcl     ~ log(1 + 0.574^2)  # IIV CL/F_ATV = 57.4% (RSE 18.1%); omega^2 = log(1 + 0.574^2)
    etalvc     ~ log(1 + 0.374^2)  # IIV V/F_ATV  = 37.4% (RSE 21.4%); omega^2 = log(1 + 0.374^2)

    etalka_rtv ~ log(1 + 2.078^2)  # IIV ka_RTV   = 207.8% (RSE 26.1%); omega^2 = log(1 + 2.078^2)
    etalcl_rtv ~ log(1 + 0.60^2)   # IIV CL/F_RTV = 60.0% (RSE 27.7%);  omega^2 = log(1 + 0.60^2)
    etalvc_rtv ~ log(1 + 0.49^2)   # IIV V/F_RTV  = 49.0% (RSE 29.1%);  omega^2 = log(1 + 0.49^2)

    # ===================================================================
    # Residual error -- ATV combined (proportional + additive); RTV
    # proportional only (the additive component of the initial combined
    # error was deleted as negligible -- Molto 2016 Results, "Basic and
    # intermediate models for RTV": "The additive component of the
    # initial combined error model used to describe the residual
    # variability was negligible, and was deleted from the model").
    # ===================================================================
    propSd     <- 0.270
    label("Atazanavir proportional residual error (fraction)")    # Table 2 final ATV residual proportional = 27.0% (RSE 7.4%)
    addSd      <- 0.07
    label("Atazanavir additive residual error (mg/L)")            # Table 2 final ATV residual additive = 0.07 mg/L (RSE 14.2%)
    propSd_rtv <- 0.280
    label("Ritonavir proportional residual error (fraction)")     # Table 2 final RTV residual proportional = 28.0% (RSE 3.5%)
  })

  model({
    # ---------------------------------------------------------------
    # 1. Ritonavir individual PK parameters. Computed first so that
    #    the RTV plasma concentration C_RTV is available for the
    #    exponential inhibition term on ATV CL/F below.
    # ---------------------------------------------------------------
    ka_rtv  <- exp(lka_rtv + etalka_rtv)
    mtt_rtv <- exp(lmtt_rtv)
    nn_rtv  <- exp(lnn_rtv)
    cl_rtv  <- exp(lcl_rtv + etalcl_rtv)
    vc_rtv  <- exp(lvc_rtv + etalvc_rtv)

    # ---------------------------------------------------------------
    # 2. Atazanavir individual PK parameters. CL/F_ATV is multiplied
    #    by an exponential RTV-concentration inhibition factor that
    #    equals 1 when C_RTV = 0 and decreases as C_RTV rises. When
    #    the RTV depot is never dosed (unboosted ATV regimen),
    #    central_rtv remains 0 and ATV CL/F equals its unboosted
    #    typical value exp(lcl) = 11.7 L/h.
    # ---------------------------------------------------------------
    crtv    <- central_rtv / vc_rtv

    ka      <- exp(lka  + etalka)
    mtt     <- exp(lmtt)
    nn      <- exp(lnn)
    cl      <- exp(lcl  + etalcl) * exp(-e_crtv_cl * crtv)
    vc      <- exp(lvc  + etalvc)

    # ---------------------------------------------------------------
    # 3. ODE system. Both drugs use rxode2's analytical Savic
    #    transit-chain (transit(n, mtt, bio)) as the input rate into
    #    the depot; transit() reads the original dose amount via
    #    podo() regardless of the f(depot) hook. The bolus
    #    contribution to depot is suppressed via f(<depot>) <- 0 so
    #    the transit() chain alone delivers the dose. Both drugs
    #    share a one-compartment-with-first-order-absorption
    #    disposition (depot -> central -> linear elimination).
    # ---------------------------------------------------------------
    d/dt(depot)           <-  transit(nn,     mtt,     1) - ka     * depot
    d/dt(central)         <-  ka     * depot     - cl     / vc     * central

    d/dt(depot_rtv)       <-  transit(nn_rtv, mtt_rtv, 1) - ka_rtv * depot_rtv
    d/dt(central_rtv)     <-  ka_rtv * depot_rtv - cl_rtv / vc_rtv * central_rtv

    # Suppress the bolus contribution; transit() provides the input.
    f(depot)     <- 0
    f(depot_rtv) <- 0

    # ---------------------------------------------------------------
    # 4. Observation variables and residual error. Cc = ATV plasma
    #    concentration (mg/L); Cc_rtv = RTV plasma concentration
    #    (mg/L). ATV uses combined additive + proportional error;
    #    RTV uses proportional error only.
    # ---------------------------------------------------------------
    Cc     <- central     / vc
    Cc_rtv <- central_rtv / vc_rtv

    Cc     ~ add(addSd) + prop(propSd)
    Cc_rtv ~ prop(propSd_rtv)
  })
}
