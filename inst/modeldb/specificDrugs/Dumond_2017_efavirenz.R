Dumond_2017_efavirenz <- function() {
  description <- "Two-compartment population PK model for efavirenz in 60 HIV-infected adults aged 22-73 years, with total and unbound (protein-free) plasma concentrations described simultaneously. The ODE system is written on the unbound side: cl, vc, q and vp are the apparent UNBOUND parameters CLu/F, Vu/F, Qu/F and Vp,u/F, the unbound concentration is Cu = central / vc, and the total concentration is Cc = Cu / fu with fu the estimated fraction unbound (0.654%). The corresponding total-drug parameters are the products CL/F = CLu/F * fu = 6.93 L/h and V/F = Vu/F * fu = 130 L. First-order absorption with no lag time. Between-occasion variability on CLu/F over two sampling occasions. Chronologic age, the Fried frailty phenotype and p16INK4a expression were screened and had no significant effect on unbound clearance or on the unbound fraction, so the final model carries no structural covariate."
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
        "$ABBR REPLACE ETA(OCC_CL)=ETA(7,8) multiplexer over the two",
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
    depot = list(analyte = "efavirenz", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "efavirenz", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "efavirenz", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 60L,
    n_studies = 1L,
    age_range = "22-73 years",
    age_median = "48 years",
    sex_female_pct = 30,
    race_ethnicity = c(`African American` = 57, White = 37, Other = 7),
    disease_state = "HIV-1 infection on stable antiretroviral therapy for at least 2 weeks; adherence at least 27 of the previous 30 doses. Median HIV duration 10.5 years (range 1-31); median CD4 count 662 cells/mm^3 (range 10-1,724). Participants with hemoglobin < 10 mg/dL, estimated creatinine clearance < 30 mL/min (Cockcroft-Gault, total body weight) or DAIDS grade 2 or higher laboratory abnormalities were excluded.",
    dose_range = "Efavirenz 600 mg by mouth once daily at steady state, co-administered with tenofovir disoproxil fumarate 300 mg and emtricitabine 200 mg once daily.",
    regions = "United States (UNC HealthCare Infectious Diseases Clinic, Chapel Hill NC; Cone Health Regional Center for Infectious Diseases, Greensboro NC). ClinicalTrials.gov NCT01180075.",
    bmi_range = "17.3-44.3 kg/m^2 (median 27.2)",
    renal_function = "Creatinine clearance (Cockcroft-Gault) median 108 mL/min (range 43-200)",
    aging_markers = "Fried frailty phenotype: 46 participants (77%) with no positive components, 8 (13%) prefrail with 1-2 components, 2 (3%) frail with 3 or more; 4 (7%) not phenotyped. Log2(p16INK4a) median 2.0 (range 0.20-2.8).",
    notes = "Table 1 (EFV arm, n = 60). 54 participants were sparsely sampled (4 samples per occasion: predose, 2 h, 4-6 h, 10-14 h postdose, over one to three occasions) and 6 provided intensive sampling (11 samples around a single observed dose). Only nonfrail participants aged 55 or older were enrolled in the intensive group. Total and unbound concentrations were measured at every sampling time (unbound by rapid equilibrium dialysis with LC-MS/MS). One EFV sample had insufficient volume for the unbound assay and was treated as missing; no total concentrations were below the limit of quantitation. NONMEM 7.3, ADVAN13, SAEM followed by importance sampling."
  )

  ini({
    # ==================================================================
    # Structural parameters - Dumond 2017 Table 2, EFV column
    # ("Estimates (RSE%)"). The paper parameterises the model on the
    # UNBOUND side: the ODE central compartment is the unbound plasma
    # compartment, so cl / vc / q / vp below are CLu/F, Vu/F, Qu/F and
    # Vp,u/F. The total-drug parameters quoted in Table 2 are the
    # products with fu (Table 2 footnotes a-d): CL/F = CLu/F * fu,
    # V/F = Vu/F * fu, Q/F = Qu/F * fu, Vp/F = Vp,u/F * fu. The
    # micro-constants are identical either way, so this is a
    # reparameterisation and not a different model.
    # Supplementary control stream PSP4-6-128-s001.txt (labelled
    # 'RTV protein binding model' in its $PROBLEM line but reading
    # efvpbnssfreeonly.csv and writing FILE=sdtabEFV, i.e. the EFV run)
    # confirms the structure: COMP(DEPOT) / COMP(FREE) / COMP(PERI),
    # CU = A(2)/V, CTOTAL = CU/FU, no ALAG1.
    # ==================================================================
    lka <- log(0.463); label("Absorption rate constant ka (1/h)") # Table 2 EFV: ka = 0.463 1/h (RSE 5%); bootstrap 0.517 (0.382-0.653)
    lcl <- log(1060); label("Apparent unbound oral clearance CLu/F (L/h)") # Table 2 EFV: CLu/F = 1,060 L/h (RSE 1%); bootstrap 1,200 (982-1,410). Total CL/F = 1,060 * 0.00654 = 6.93 L/h (Table 2 row 'CL/F', footnote a)
    lvc <- log(19900); label("Apparent unbound central volume Vu/F (L)") # Table 2 EFV: Vu/F = 19,900 L (RSE 1%); bootstrap 20,700 (14,800-29,700). Total V/F = 19,900 * 0.00654 = 130 L (Table 2 row 'V/F', footnote b)
    lq <- log(5940); label("Apparent unbound intercompartmental clearance Qu/F (L/h)") # Table 2 EFV: Qu/F = 5,940 L/h (RSE 2%); bootstrap 3,980 (2,370-7,940). Total Q/F = 5,940 * 0.00654 = 38.8 L/h (Table 2 row 'Q/F', footnote c)
    lvp <- log(24300); label("Apparent unbound peripheral volume Vp,u/F (L)") # Table 2 EFV: Vp,u/F = 24,300 L (RSE 1%); bootstrap 22,000 (18,600-28,400). Total Vp/F = 24,300 * 0.00654 = 159 L (Table 2 row 'Vp/F', footnote d)
    lfu <- log(0.00654); label("Fraction of efavirenz unbound in plasma fu (unitless)") # Table 2 EFV: fu = 0.654% (RSE 0.5%); bootstrap 0.647% (0.621-0.681%). Entered as the fraction 0.00654

    # ==================================================================
    # Inter-individual variability - Dumond 2017 Table 2, first
    # 'IIV (CV%)' block of the EFV column. Methods: IIV was 'assumed to
    # be normally distributed and exponentially related to the
    # population parameters', i.e. P_i = TVP * exp(eta_i), so the
    # log-scale variance is omega^2 = log(1 + CV^2). Shrinkage as
    # printed in brackets by Table 2.
    #
    # DEVIATION: the supplementary control stream fits a full
    # $OMEGA BLOCK(6) across all six etas. Table 2 prints only the
    # diagonal CV% values, so the off-diagonal covariances are not
    # recoverable and the IIV is encoded here as diagonal. See the
    # vignette's Assumptions and deviations section.
    # ==================================================================
    etalcl ~ 0.0812859 # Table 2 EFV IIV CLu/F = 29.1% CV (RSE 1%) [shrinkage 19%]; bootstrap 33.7 (6.91-51.9) -> log(1 + 0.291^2)
    etalvc ~ 0.385668 # Table 2 EFV IIV Vu/F = 68.6% CV (RSE 1%) [shrinkage 19%]; bootstrap 47.1 (3.98-103) -> log(1 + 0.686^2)
    etalka ~ 0.0249667 # Table 2 EFV IIV ka = 15.9% CV (RSE 3%) [shrinkage 18%]; bootstrap 19.0 (3.12-67.4) -> log(1 + 0.159^2)
    etalfu ~ 0.0301654 # Table 2 EFV IIV fu = 17.5% CV (RSE 14%) [shrinkage 14%]; bootstrap 16.2 (9.54-22.4) -> log(1 + 0.175^2)
    etalq ~ 0.218361 # Table 2 EFV IIV Qu/F = 49.4% CV (RSE 8%) [shrinkage 21%]; bootstrap 43.8 (6.66-143) -> log(1 + 0.494^2)
    etalvp ~ 0.0969098 # Table 2 EFV IIV Vp,u/F = 31.9% CV (RSE 2%) [shrinkage 19%]; bootstrap 17.7 (3.47-45.7) -> log(1 + 0.319^2)

    # ==================================================================
    # Between-occasion variability on CLu/F. Results: 'Interoccasional
    # variability was incorporated on CLu/F for the three drugs.'
    # Table 2 prints this as a SECOND block also headed 'IIV (CV%)'
    # containing only a CLu/F row; that heading is a typesetting error
    # in the published table - the block is the interoccasion variance,
    # matching the control stream's $OMEGA BLOCK(1) + BLOCK(1) SAME
    # pair multiplexed by $ABBR REPLACE ETA(OCC_CL)=ETA(7,8).
    # The SAME keyword means both occasions share one variance, so the
    # second occasion is fixed to the first estimate.
    # ==================================================================
    etaiov_cl_1 ~ 0.165325 # Table 2 EFV second 'IIV (CV%)' block, CLu/F = 42.4% CV (RSE 36%); bootstrap 44.0 (22.7-60.1) -> log(1 + 0.424^2)
    etaiov_cl_2 ~ fix(0.165325) # control stream $OMEGA BLOCK(1) SAME: occasion 2 shares the occasion-1 variance

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
    propSd <- 0.242; label("Proportional residual error on total plasma efavirenz Cc (fraction)") # Table 2 EFV residual error, Total = 24.2% CV (RSE 10%); bootstrap 24.9 (19.8-30.8)
    propSd_Cu <- 0.313; label("Proportional residual error on unbound plasma efavirenz Cu (fraction)") # Table 2 EFV residual error, Unbound = 31.3% CV (RSE 20%); bootstrap 30.9 (25.1-37.0)
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
    # FU = EXP(MU_4 + ETA(4)), Q = EXP(MU_5 + ETA(5)),
    # VP = EXP(MU_6 + ETA(6)). No covariate enters the final EFV model.
    cl <- exp(lcl + etalcl + iov_cl)
    vc <- exp(lvc + etalvc)
    ka <- exp(lka + etalka)
    fu <- exp(lfu + etalfu)
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

    # ---- Residual error ---------------------------------------------
    Cc ~ prop(propSd)
    Cu ~ prop(propSd_Cu)
  })
}
