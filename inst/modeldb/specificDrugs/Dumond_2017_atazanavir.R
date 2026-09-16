Dumond_2017_atazanavir <- function() {
  description <- "Two-compartment population PK model with a first-order absorption lag time for ritonavir-boosted atazanavir in 31 HIV-infected adults aged 24-61 years, with total and unbound (protein-free) plasma concentrations described simultaneously. The ODE system is written on the unbound side: cl, vc, q and vp are the apparent UNBOUND parameters CLu/F, Vu/F, Qu/F and Vp,u/F, the unbound concentration is Cu = central / vc, and the total concentration is Cc = Cu / fu with fu the estimated fraction unbound (5.67%). The corresponding total-drug parameters are the products CL/F = CLu/F * fu = 5.73 L/h and V/F = Vu/F * fu = 62.9 L. Between-occasion variability on CLu/F over two sampling occasions. Chronologic age, the Fried frailty phenotype, p16INK4a expression, body weight and body mass index were screened and none was retained for atazanavir, so the final model carries no structural covariate. The companion ritonavir model from the same paper is Dumond_2017_ritonavir."
  reference <- paste(
    "Dumond JB, Chen J, Cottrell M, Trezza CR, Prince HMA, Sykes C, Torrice C,",
    "White N, Malone S, Wang R, Patterson KB, Sharpless NE, Forrest A.",
    "Population pharmacokinetics modeling of unbound efavirenz, atazanavir,",
    "and ritonavir in HIV-infected subjects with aging biomarkers.",
    "CPT Pharmacometrics Syst Pharmacol. 2017;6(2):128-135.",
    "doi:10.1002/psp4.12151.",
    sep = " "
  )
  vignette <- "Dumond_2017_unbound_antiretrovirals"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    OCC = list(
      description = "Sampling-occasion index used for the between-occasion random effect on unbound clearance",
      units = "(index)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Two occasions, matching the supplementary control stream's",
        "$ABBR REPLACE ETA(OCC_CL)=ETA(8,9) multiplexer over the two",
        "$OMEGA BLOCK(1) / BLOCK(1) SAME interoccasion etas on CLu/F.",
        "Sparse-sampling participants contributed one to three occasions in",
        "the study; the final control stream allocates two IOV etas. For a",
        "single-occasion simulation set OCC = 1 throughout.",
        sep = " "
      ),
      source_name = "OCC"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "atazanavir", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "atazanavir", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "atazanavir", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 31L,
    n_studies = 1L,
    age_range = "24-61 years",
    age_median = "49 years",
    sex_female_pct = 39,
    race_ethnicity = c(`African American` = 61, White = 32, Other = 6),
    disease_state = "HIV-1 infection on stable antiretroviral therapy for at least 2 weeks; adherence at least 27 of the previous 30 doses. Median HIV duration 10 years (range 1-24); median CD4 count 692 cells/mm^3 (range 375-1,501). Participants with hemoglobin < 10 mg/dL, estimated creatinine clearance < 30 mL/min (Cockcroft-Gault, total body weight) or DAIDS grade 2 or higher laboratory abnormalities were excluded, except that drug-attributed total bilirubin elevations were permitted.",
    dose_range = "Atazanavir 300 mg with ritonavir 100 mg by mouth once daily at steady state, co-administered with tenofovir disoproxil fumarate 300 mg and emtricitabine 200 mg once daily.",
    regions = "United States (UNC HealthCare Infectious Diseases Clinic, Chapel Hill NC; Cone Health Regional Center for Infectious Diseases, Greensboro NC). ClinicalTrials.gov NCT01180075.",
    bmi_range = "20.2-40.4 kg/m^2 (median 30.3)",
    renal_function = "Creatinine clearance (Cockcroft-Gault) median 100 mL/min (range 67-227)",
    aging_markers = "Fried frailty phenotype: 21 participants (68%) with no positive components, 9 (29%) prefrail with 1-2 components, 1 (3%) frail with 3 or more. Log2(p16INK4a) median 2.2 (range 0.16-3.9).",
    notes = "Table 1 (ATV/RTV arm, n = 31). 25 participants were sparsely sampled (4 samples per occasion: predose, 2 h, 4-6 h, 10-14 h postdose, over one to three occasions) and 6 provided intensive sampling (11 samples around a single observed dose). Only nonfrail participants aged 55 or older were enrolled in the intensive group. Total and unbound concentrations were measured at every sampling time (unbound by rapid equilibrium dialysis with LC-MS/MS). One ATV sample was below the limit of quantitation and was imputed as half the lower limit. NONMEM 7.3, ADVAN13, SAEM followed by importance sampling. Because of sparse sampling and substantial absorption-phase variability, ka and the lag time were estimated with RSE > 30%."
  )

  ini({
    # ==================================================================
    # Structural parameters - Dumond 2017 Table 2, ATV column
    # ("Estimates (RSE%)"). The paper parameterises the model on the
    # UNBOUND side: the ODE central compartment is the unbound plasma
    # compartment, so cl / vc / q / vp below are CLu/F, Vu/F, Qu/F and
    # Vp,u/F. The total-drug parameters quoted in Table 2 are the
    # products with fu (Table 2 footnotes a-d): CL/F = CLu/F * fu,
    # V/F = Vu/F * fu, Q/F = Qu/F * fu, Vp/F = Vp,u/F * fu. The
    # micro-constants are identical either way, so this is a
    # reparameterisation and not a different model.
    # Supplementary control stream PSP4-6-128-s002.txt
    # ('ATV protein binding model with free as core') confirms the
    # structure: COMP(DEPOT) / COMP(FREE) / COMP(PERI), ALAG1 on the
    # depot, CU = A(2)/V, CTOTAL = CU/FU.
    # ==================================================================
    lka <- log(0.705); label("Absorption rate constant ka (1/h)") # Table 2 ATV: ka = 0.705 1/h (RSE 42%); bootstrap 0.623 (0.452-0.885)
    ltlag <- log(0.529); label("Absorption lag time (h)") # Table 2 ATV: lag time = 0.529 h (RSE 36%); bootstrap 0.577 (0.165-1.10)
    lcl <- log(101); label("Apparent unbound oral clearance CLu/F (L/h)") # Table 2 ATV: CLu/F = 101 L/h (RSE 2%); bootstrap 105 (79.0-142). Total CL/F = 101 * 0.0567 = 5.73 L/h (Table 2 row 'CL/F', footnote a)
    lvc <- log(1110); label("Apparent unbound central volume Vu/F (L)") # Table 2 ATV: Vu/F = 1,110 L (RSE 3%); bootstrap 1,090 (778-1,620). Total V/F = 1,110 * 0.0567 = 62.9 L (Table 2 row 'V/F', footnote b)
    lq <- log(134); label("Apparent unbound intercompartmental clearance Qu/F (L/h)") # Table 2 ATV: Qu/F = 134 L/h (RSE 1%); bootstrap 137 (91.4-183). Total Q/F = 134 * 0.0567 = 7.61 L/h (Table 2 row 'Q/F', footnote c)
    lvp <- log(2280); label("Apparent unbound peripheral volume Vp,u/F (L)") # Table 2 ATV: Vp,u/F = 2,280 L (RSE 1%); bootstrap 2,440 (1,730-2,950). Total Vp/F = 2,280 * 0.0567 = 129 L (Table 2 row 'Vp/F', footnote d)
    lfu <- log(0.0567); label("Fraction of atazanavir unbound in plasma fu (unitless)") # Table 2 ATV: fu = 5.67% (RSE 1%); bootstrap 5.73% (5.34-6.11%). Entered as the fraction 0.0567

    # ==================================================================
    # Inter-individual variability - Dumond 2017 Table 2, first
    # 'IIV (CV%)' block of the ATV column. Methods: IIV was 'assumed to
    # be normally distributed and exponentially related to the
    # population parameters', i.e. P_i = TVP * exp(eta_i), so the
    # log-scale variance is omega^2 = log(1 + CV^2). Shrinkage as
    # printed in brackets by Table 2.
    #
    # DEVIATION: the supplementary control stream fits a full
    # $OMEGA BLOCK(7) across all seven etas. Table 2 prints only the
    # diagonal CV% values, so the off-diagonal covariances are not
    # recoverable and the IIV is encoded here as diagonal. See the
    # vignette's Assumptions and deviations section.
    # ==================================================================
    etalcl ~ 0.128305 # Table 2 ATV IIV CLu/F = 37.0% CV (RSE 1%) [shrinkage 8%]; bootstrap 35.1 (26.4-44.0) -> log(1 + 0.370^2)
    etalvc ~ 0.179191 # Table 2 ATV IIV Vu/F = 44.3% CV (RSE 1%) [shrinkage 14%]; bootstrap 42.3 (31.0-72.0) -> log(1 + 0.443^2)
    etalka ~ 0.113075 # Table 2 ATV IIV ka = 34.6% CV (RSE 1%) [shrinkage 11%]; bootstrap 35.0 (20.0-54.5) -> log(1 + 0.346^2)
    etaltlag ~ 0.500555 # Table 2 ATV IIV lag time = 80.6% CV (RSE 8%) [shrinkage 19%]; bootstrap 77.0 (58.6-104) -> log(1 + 0.806^2)
    etalfu ~ 0.0180609 # Table 2 ATV IIV fu = 13.5% CV (RSE 4%) [shrinkage 7%]; bootstrap 14.2 (9.1-18.5) -> log(1 + 0.135^2)
    etalq ~ 0.0502453 # Table 2 ATV IIV Qu/F = 22.7% CV (RSE 14%) [shrinkage 23%]; bootstrap 20.1 (3.82-54.6) -> log(1 + 0.227^2)
    etalvp ~ 0.0498144 # Table 2 ATV IIV Vp,u/F = 22.6% CV (RSE 13.0%) [shrinkage 30%]; bootstrap 19.8 (3.84-54.0) -> log(1 + 0.226^2)

    # ==================================================================
    # Between-occasion variability on CLu/F. Results: 'Interoccasional
    # variability was incorporated on CLu/F for the three drugs.'
    # Table 2 prints this as a SECOND block also headed 'IIV (CV%)'
    # containing only a CLu/F row; that heading is a typesetting error
    # in the published table - the block is the interoccasion variance,
    # matching the control stream's $OMEGA BLOCK(1) + BLOCK(1) SAME
    # pair multiplexed by $ABBR REPLACE ETA(OCC_CL)=ETA(8,9).
    # The SAME keyword means both occasions share one variance, so the
    # second occasion is fixed to the first estimate.
    # ==================================================================
    etaiov_cl_1 ~ 0.435749 # Table 2 ATV second 'IIV (CV%)' block, CLu/F = 73.9% CV (RSE 66%); bootstrap 69.2 (39.4-100) -> log(1 + 0.739^2)
    etaiov_cl_2 ~ fix(0.435749) # control stream $OMEGA BLOCK(1) SAME: occasion 2 shares the occasion-1 variance

    # ==================================================================
    # Residual error - Dumond 2017 Table 2 'Residual error (CV%)'.
    # Control stream $ERROR: W = IPRED and
    # Y = IPRED + W*EPS(1)*(1-TYPE) + W*EPS(2)*TYPE with TYPE = 0 for
    # DVID 1 (unbound) and TYPE = 1 for DVID 2 (total), i.e. a
    # proportional error with a separate magnitude per analyte stream.
    #
    # DEVIATION: the control stream's $SIGMA BLOCK(2) estimates the
    # correlation between the two residual streams (the paper's L2 data
    # item). Only the two diagonal CV% values are published, and
    # nlmixr2 has no cross-endpoint residual correlation, so the two
    # streams are encoded as independent here.
    # ==================================================================
    propSd <- 0.276; label("Proportional residual error on total plasma atazanavir Cc (fraction)") # Table 2 ATV residual error, Total = 27.6% CV (RSE 10%); bootstrap 27.8 (23.8-32.4)
    propSd_Cu <- 0.301; label("Proportional residual error on unbound plasma atazanavir Cu (fraction)") # Table 2 ATV residual error, Unbound = 30.1% CV (RSE 11%); bootstrap 30.6 (26.9-34.8)
  })

  model({
    # ---- Between-occasion variability multiplexer -------------------
    # Control stream: ETA(OCC_CL) selects one of two interoccasion etas
    # by the OCC data item; the selected eta is added to the CLu/F
    # log-scale typical value alongside the subject-level eta.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2

    # ---- Individual parameters (unbound-side parameterisation) ------
    # Control stream $PK: CL = EXP(MU_1 + ETA(1) + ETA(OCC_CL)),
    # V = EXP(MU_2 + ETA(2)), KA = EXP(MU_3 + ETA(3)),
    # FU = EXP(MU_4 + ETA(4)), ALAG1 = EXP(MU_5 + ETA(5)),
    # Q = EXP(MU_6 + ETA(6)), VP = EXP(MU_7 + ETA(7)).
    # No covariate enters the final ATV model.
    cl <- exp(lcl + etalcl + iov_cl)
    vc <- exp(lvc + etalvc)
    ka <- exp(lka + etalka)
    fu <- exp(lfu + etalfu)
    tlag <- exp(ltlag + etaltlag)
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)

    # ---- Concentrations --------------------------------------------
    # Control stream $ERROR: CU = A(2)/V and CTOTAL = CU/FU. The
    # central state holds the amount of drug; dividing by the apparent
    # unbound volume gives the unbound concentration and dividing that
    # by fu gives the total concentration. Equivalently
    # Cc = central / (vc * fu) = central / (V/F).
    Cu <- central / vc
    Cc <- Cu / fu

    # ---- ODE system -------------------------------------------------
    # Control stream $DES:
    #   DADT(1) = -KA*A(1)
    #   DADT(2) =  KA*A(1) - CL/V*A(2) - Q/V*A(2) + Q/VP*A(3)
    #   DADT(3) =  Q/V*A(2) - Q/VP*A(3)
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - cl * Cu - q * (Cu - peripheral1 / vp)
    d/dt(peripheral1) <- q * (Cu - peripheral1 / vp)

    # ---- Absorption lag ---------------------------------------------
    # Control stream ALAG1 targets compartment 1 = COMP(DEPOT).
    alag(depot) <- tlag

    # ---- Residual error ---------------------------------------------
    Cc ~ prop(propSd)
    Cu ~ prop(propSd_Cu)
  })
}
