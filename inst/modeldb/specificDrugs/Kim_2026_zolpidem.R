Kim_2026_zolpidem <- function() {
  description <- paste(
    "One-compartment population PK model for oral immediate-release",
    "zolpidem in 30 healthy Korean volunteers (15 male, 15 female) given a",
    "single 10 mg dose (Kim 2026). Absorption is a Savic transit-compartment",
    "chain feeding the `depot` state, which then empties first-order at ka",
    "into `central`; disposition is one-compartment with first-order",
    "elimination. The transit chain is parameterised by the mean transit",
    "time (MTT 0.25 h) and the number of transit compartments (NN 19.4),",
    "with ktr = (NN + 1) / MTT. Relative bioavailability is fixed at a",
    "typical value of 1 but carries estimated inter-individual variability;",
    "the authors added it specifically to absorb the CL/F-Vc/F covariance,",
    "after which the IIV on Vc/F estimated to zero. No covariate met the",
    "retention criteria in the final PK model: sex entered on forward",
    "selection but was dropped on backward elimination. This is the PK",
    "backbone that the three companion PK-PD models",
    "(`Kim_2026_zolpidem_dsst`, `Kim_2026_zolpidem_crt`,",
    "`Kim_2026_zolpidem_vas`) each hold fixed."
  )
  reference <- paste(
    "Kim HC, Sunwoo J, Yoon S, Jang I-J, Chung J-Y. Population",
    "Pharmacokinetic-Pharmacodynamic Analysis of Zolpidem in Healthy",
    "Volunteers: Covariate Contributions to Variability in CNS Responses.",
    "CPT Pharmacometrics Syst Pharmacol. 2026;15:e70208.",
    "doi:10.1002/psp4.70208"
  )
  vignette <- "Kim_2026_zolpidem"
  units <- list(time = "h", dosing = "mg", concentration = "ug/L")

  compartmentData <- list(
    depot   = list(analyte = "zolpidem", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "zolpidem", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list()

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Sex, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "Kim 2026 Results 3.2.1: sex was 'the only covariate selected",
        "during the forward selection step (dOFV = -5.137) but was removed",
        "during the subsequent backward elimination step and was therefore",
        "not retained in the final PK model' (Tables S3 and S4). The",
        "Discussion attributes this to limited power rather than to absence",
        "of a true effect: the companion NCA reported a 22% lower CL/F in",
        "females that was not statistically significant. No point estimate",
        "is reported for the forward-selection step, so the effect cannot",
        "be encoded. Source column SEX coded 1 = M, 2 = F."
      ),
      source_name = "SEX"
    ),
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      notes = paste(
        "Screened in the stepwise covariate analysis of the PK model but",
        "not retained (Kim 2026 Results 3.2.1, Tables S3 and S4). Cohort",
        "median 60.05 kg (range 50.0-83.0). Weight IS retained on the PD",
        "potency parameters of the companion DSST and CRT models. The",
        "Discussion notes that adding weight on the IIV of bioavailability",
        "reduced the OFV by 4.554, below the forward-selection threshold.",
        "Source column WEI."
      ),
      source_name = "WEI"
    ),
    AGE = list(
      description = "Age at baseline",
      units = "years",
      type = "continuous",
      notes = paste(
        "Screened but not retained in the final PK model (Kim 2026 Results",
        "3.2.1). Cohort median 29.5 years (range 20-44). Age IS retained on",
        "the PD baselines of the companion DSST and CRT models."
      ),
      source_name = "AGE"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/dL as reported by Kim 2026 (canonical register unit is g/L)",
      type = "continuous",
      notes = paste(
        "Screened but not retained in the final PK model (Kim 2026 Results",
        "3.2.1). Cohort median 4.4 g/dL (range 3.9-5.0). Albumin IS",
        "retained on the PD potency and Hill parameters of the companion",
        "DSST and CRT models."
      ),
      source_name = "ALB"
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "IU/L",
      type = "continuous",
      notes = paste(
        "Screened in the stepwise covariate analysis but not retained in",
        "any final model (Kim 2026 Results 3.2.1, Table S1). Cohort median",
        "13.5 IU/L (range 6-46). No point estimate reported."
      ),
      source_name = "ALT"
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "IU/L",
      type = "continuous",
      notes = paste(
        "Screened in the stepwise covariate analysis but not retained in",
        "any final model (Kim 2026 Results 3.2.1, Table S1). Cohort median",
        "22 IU/L (range 12-42). No point estimate reported."
      ),
      source_name = "AST"
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "mg/dL as reported by Kim 2026 (canonical register unit is umol/L)",
      type = "continuous",
      notes = paste(
        "Screened in the stepwise covariate analysis but not retained in",
        "the final PK model (Kim 2026 Results 3.2.1, Table S1). Cohort",
        "median 0.655 mg/dL (range 0.32-1.42). Total bilirubin reached",
        "statistical significance on H1 of the companion VAS model but was",
        "excluded there as biologically implausible (Results 3.2.4). No",
        "point estimate reported. Source column TBIL."
      ),
      source_name = "TBIL"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      notes = paste(
        "Screened in the stepwise covariate analysis but not retained in",
        "any final model (Kim 2026 Results 3.2.1, Table S1). Cohort median",
        "0.73 mg/dL (range 0.45-0.97). No point estimate reported. Source",
        "column CREA."
      ),
      source_name = "CREA"
    ),
    DSST_BL = list(
      description = "Baseline digit symbol substitution test score at 1 h on Day -1",
      units = "score",
      type = "continuous",
      notes = paste(
        "Screened as a candidate covariate on the PK parameters but not",
        "retained (Kim 2026 Methods 2.4, Table S1; control-stream $INPUT",
        "column DSSTB1, 'DSST baseline 1h'). Cohort median 72.5 (range",
        "46-104). No point estimate reported."
      ),
      source_name = "DSSTB1"
    ),
    CRT_BL = list(
      description = "Baseline choice reaction time at 1 h on Day -1",
      units = "msec",
      type = "continuous",
      notes = paste(
        "Screened as a candidate covariate on the PK parameters but not",
        "retained (Kim 2026 Methods 2.4, Table S1; control-stream $INPUT",
        "column CRTB1, 'CRT baseline 1h'). Cohort median 438.5 msec (range",
        "372-709). CRT_BL IS retained on EC50 and HILL in the companion",
        "PK-CRT model."
      ),
      source_name = "CRTB1"
    ),
    VAS_SEDATION_BL = list(
      description = "Baseline sedation visual analog scale at 1 h on Day -1",
      units = "mm",
      type = "continuous",
      notes = paste(
        "Screened as a candidate covariate on the PK parameters but not",
        "retained (Kim 2026 Methods 2.4, Table S1; control-stream $INPUT",
        "column VASB1, 'VAS baseline 1h'). Cohort median 30.5 mm (range",
        "0-97). No point estimate reported."
      ),
      source_name = "VASB1"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 30,
    n_studies = 1,
    age_range = "20-44 years",
    age_median = "29.5 years",
    weight_range = "50.0-83.0 kg",
    weight_median = "60.05 kg",
    sex_female_pct = 50,
    race_ethnicity = c(Asian = 100),
    disease_state = "healthy volunteers",
    dose_range = "single oral 10 mg immediate-release zolpidem",
    regions = "Korea",
    n_observations = "325 plasma zolpidem concentrations",
    laboratory = paste(
      "Median (range) at baseline (Kim 2026 Table 1): albumin 4.4",
      "(3.9-5.0) g/dL; ALT 13.5 (6-46) IU/L; AST 22 (12-42) IU/L; total",
      "bilirubin 0.655 (0.32-1.42) mg/dL; serum creatinine 0.73",
      "(0.45-0.97) mg/dL."
    ),
    notes = paste(
      "Secondary population PK-PD analysis of a previously reported",
      "open-label single-dose study in healthy Korean volunteers. PK",
      "sampling at 0, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8 and 12 h on",
      "Day 1 (Kim 2026 Methods 2.1). Baseline characteristics are given in",
      "Kim 2026 Table 1; weight, albumin, ALT, AST, total bilirubin and",
      "serum creatinine all differed significantly between sexes.",
      "Concentrations below the limit of quantification were handled by",
      "Beal's M1 method (records set MDV = 1 and dropped)."
    )
  )

  ini({
    # ====================================================================
    # STRUCTURAL PHARMACOKINETIC PARAMETERS
    # ====================================================================
    # Typical values are the final-model estimates of Kim 2026 Table 2
    # ("Population pharmacokinetic parameter estimates and sampling
    # importance resampling results"), Final model Estimates (%RSE)
    # column, and are reproduced verbatim in the $THETA block of Data S1
    # (the deposited NONMEM control stream for the final PK model).
    # CL/F and Vc/F are apparent (oral) parameters: bioavailability is
    # fixed at 1 rather than estimated.

    lka <- log(11.7)
    label("Absorption rate constant, depot -> central (Ka, 1/h)")                    # Kim 2026 Table 2: Ka = 11.7 1/h (RSE 36.8%; SIR median 12.1, 95% CI 6.6-22.3). Data S1 $THETA(3).

    lcl <- log(18.0)
    label("Apparent clearance (CL/F, L/h)")                                          # Kim 2026 Table 2: CL/F = 18.0 L/h (RSE 7.9%; SIR median 18.0, 95% CI 15.9-20.2). Data S1 $THETA(1).

    lvc <- log(64.0)
    label("Apparent central volume of distribution (Vc/F, L)")                       # Kim 2026 Table 2: Vc/F = 64.0 L (RSE 4.9%; SIR median 64.1, 95% CI 58.3-69.9). Data S1 $THETA(2).

    lmtt <- log(0.25)
    label("Mean transit time of the absorption chain (MTT, h)")                      # Kim 2026 Table 2: MTT = 0.25 h (RSE 8.4%; SIR median 0.25, 95% CI 0.21-0.29). Data S1 $THETA(4).

    lnn <- log(19.4)
    label("Number of transit compartments (NN, unitless)")                           # Kim 2026 Table 2: NN = 19.4 (RSE 18.8%; SIR median 19.4, 95% CI 13.6-28.9). Data S1 $THETA(5). Estimated as a continuous quantity, so the absorption input is the closed-form Erlang density rather than an integer chain of ODE states.

    lfdepot <- fixed(log(1))
    label("Relative bioavailability (BIO, fraction)")                                # Kim 2026 Table 2: BIO = "1 Fixed". Data S1 $THETA(6) is "(1) FIX". Results 3.2.1: bioavailability "was included in the model with a fixed typical value of 1 and estimated IIV, resulting in an improvement in the OFV (dOFV = -29.938) and a reduction in the IIV of apparent clearance (CL/F) and apparent volume of distribution (Vc/F)".

    # ====================================================================
    # INTER-INDIVIDUAL VARIABILITY
    # ====================================================================
    # Data S1 $OMEGA is a diagonal block of VARIANCES. Kim 2026 Table 2
    # reports them as %CV; the note under Table 3 gives the transform
    # used, %CV = 100 * sqrt(exp(omega^2) - 1), and every reported %CV
    # reproduces exactly from the control-stream variance under it. IIV on
    # Vc/F and on NN were both estimated to zero and fixed there, so no
    # etalvc / etalnn term appears here.

    etalka ~ 2.48
    # Kim 2026 Table 2: IIV Ka = 330.8 %CV (RSE 30.7%, shrinkage 27.9%). Data S1 $OMEGA(3,3) = 2.48; 100*sqrt(exp(2.48)-1) = 330.8%. Results 3.2.1 flags this as "relatively high", reflecting rapid and variable absorption of a BCS class I drug.

    etalcl ~ 0.0565
    # Kim 2026 Table 2: IIV CL/F = 24.1 %CV (RSE 28.1%, shrinkage 4.70%). Data S1 $OMEGA(1,1) = 0.0565; 100*sqrt(exp(0.0565)-1) = 24.1%.

    etalmtt ~ 0.133
    # Kim 2026 Table 2: IIV MTT = 37.7 %CV (RSE 51.0%, shrinkage 20.0%). Data S1 $OMEGA(4,4) = 0.133; 100*sqrt(exp(0.133)-1) = 37.7%.

    etalfdepot ~ 0.0647
    # Kim 2026 Table 2: IIV BIO = 25.9 %CV (RSE 25.3%, shrinkage 5.50%). Data S1 $OMEGA(6,6) = 0.0647; 100*sqrt(exp(0.0647)-1) = 25.9%. The typical value is fixed at 1 but the IIV is estimated -- the term exists to carry the CL/F-Vc/F covariance.

    # ====================================================================
    # RESIDUAL UNEXPLAINED VARIABILITY
    # ====================================================================
    # Data S1 $ERROR builds a combined variance,
    #   W = SQRT(THETA(7)**2 * IPRED**2 + THETA(8)**2),
    # but THETA(8) is fixed at 1e-5 ug/L, i.e. a numerical guard against a
    # zero standard deviation rather than an estimated additive
    # component. Kim 2026 Table 2 accordingly reports only the
    # proportional term and Results 3.2.1 describes the model as
    # "a proportional error model".

    propSd <- 0.189
    label("Proportional residual error (fraction)")                                  # Kim 2026 Table 2: proportional residual error = 18.9% (RSE 7.20%; SIR median 18.9, 95% CI 17.3-20.8). Data S1 $THETA(7) = 0.189.

    addSd <- fixed(0.00001)
    label("Additive residual error (ug/L)")                                          # Data S1 $THETA(8) = "(0, 0.00001) FIX". Not reported in Kim 2026 Table 2; retained here so the encoded error model matches the deposited $ERROR block exactly.
  })

  model({
    # --- 1. Individual parameters ---------------------------------------
    # No covariate was retained in the final PK model, so each individual
    # parameter is the typical value times its exponential IIV term.
    ka     <- exp(lka + etalka)
    cl     <- exp(lcl + etalcl)
    vc     <- exp(lvc)
    mtt    <- exp(lmtt + etalmtt)
    nn     <- exp(lnn)
    fdepot <- exp(lfdepot + etalfdepot)

    # --- 2. Micro-constants ----------------------------------------------
    kel <- cl / vc

    # --- 3. ODE system ----------------------------------------------------
    # Data S1 declares COMP(DEPOT) and COMP(CENTRAL) and integrates
    #   DADT(1) = INP - KA * A(1)
    #   DADT(2) = KA * A(1) - K * A(2)
    # where INP is the Savic 2007 closed-form transit input written out
    # in $DES as
    #   KTR = (NN + 1) / MTT
    #   LNFAC = LOG(2.5066) + (NN + 0.5) * LOG(NN) - NN
    #   INP = EXP(LOG(BIO*PODO) + LOG(KTR) + NN*LOG(KTR*TAD)
    #             - KTR*TAD - LNFAC)
    # rxode2's builtin transit(nn, mtt, fdepot) evaluates exactly this
    # density, deriving ktr internally as (nn + 1) / mtt and normalising
    # with the exact lgamma(nn + 1) in place of the authors' Stirling
    # approximation LNFAC (the two differ by ~0.4% of the input rate at
    # NN = 19.4). f(depot) <- 0 below suppresses the ordinary dose bolus so
    # the transit density is the only input pathway.
    d/dt(depot)   <- transit(nn, mtt, fdepot) - ka * depot
    d/dt(central) <- ka * depot - kel * central

    # --- 4. Bioavailability ----------------------------------------------
    # Data S1 $PK sets F1 = 0 for the same reason.
    f(depot) <- 0

    # --- 5. Observation and residual error --------------------------------
    # Data S1 scales the central amount with S2 = V2/1000, i.e. a dose in
    # mg over a volume in L reported as ug/L.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd) + add(addSd)
  })
}
