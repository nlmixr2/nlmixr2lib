Hodiamont_2017_gentamicin <- function() {
  description <- "Two-compartment population PK model of intravenous gentamicin in critically ill adult ICU patients (Hodiamont 2017) estimated without retained covariates, with correlated between-subject variability on CL and central volume V1, combined additive plus proportional residual error, and substantial inter-occasion variability on CL and V1 reported in the source (documented in the vignette assumptions, not encoded structurally)."
  reference <- paste(
    "Hodiamont CJ, Janssen JM, de Jong MD, Mathot RA, Juffermans NP, van Hest RM.",
    "Therapeutic Drug Monitoring of Gentamicin Peak Concentrations in Critically Ill Patients.",
    "Ther Drug Monit 2017;39(5):522-530.",
    "doi:10.1097/FTD.0000000000000432.",
    sep = " "
  )
  vignette <- "Hodiamont_2017_gentamicin"
  units    <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  # Hodiamont 2017 tested total, ideal, and adjusted body weight as covariates
  # on CL, Q, V1, and V2 using both allometric and univariate forms, and
  # reported that no body-weight covariate improved the OFV by 3.84 or more
  # (Results "Model" paragraph). Renal function, fluid balance, and inotrope
  # use were intentionally not tested; the authors wanted to quantify the
  # total inter-occasion variability seen in routine TDM where only body
  # weight is normally considered. The final model therefore retains no
  # covariate effects, so covariateData is empty.
  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 59L,
    n_studies      = 1L,
    age_mean       = "60.9 +/- 17.2 years",
    weight_mean    = "TBW 79.2 +/- 22.0 kg; IBW 71.4 +/- 11.6 kg; ABW 74.6 +/- 13.0 kg",
    sex_female_pct = 49,
    disease_state  = paste(
      "Critically ill adults admitted to a mixed medical-surgical ICU",
      "receiving intravenous gentamicin. Four of 62 treatment episodes",
      "were for endocarditis (3 mg/kg combined with a beta-lactam for",
      "synergy); these were included in the PK model fit but excluded",
      "from the primary TDM end points."
    ),
    dose_range     = paste(
      "Fixed first dose approximately 5 mg/kg (mean 5.1 +/- 1.1 mg/kg",
      "TBW for the 58 non-endocarditis episodes), 3 mg/kg for the 4",
      "endocarditis episodes, administered as a 30-min IV infusion.",
      "Mean of 2.1 +/- 1.9 administrations per treatment episode."
    ),
    n_observations = 416L,
    n_administrations = 130L,
    n_treatment_episodes = 62L,
    crcl_summary = paste(
      "Cockcroft-Gault CL: 87.0 +/- 64.7 mL/min before the first dose",
      "(n = 62 episodes), 99.7 +/- 59.3 before the second dose (n = 33),",
      "133.0 +/- 85.6 before the third (n = 13)."
    ),
    regions        = "Single-centre cohort: Academic Medical Center ICU, Amsterdam, the Netherlands. Data collected May-June 2013 and April-June 2014.",
    notes          = paste(
      "Baseline demographics from Hodiamont 2017 Table 1. Renal function",
      "showed large between-patient and within-patient variation; in part",
      "this is what the unmodelled IOV on CL captures. The cohort was",
      "fitted in NONMEM 7.2 with FOCE-I."
    )
  )

  ini({
    # Structural population PK parameters (Hodiamont 2017 Table 2, Final
    # Model column). All four typical values are reported on the linear
    # scale; nlmixr2lib log-transforms positive-constrained parameters.
    lcl <- log(2.3);   label("Clearance (CL, L/h)")                                  # Table 2: CL = 2.3 L/h
    lvc <- log(21.6);  label("Central volume of distribution (V1, L)")               # Table 2: V1 = 21.6 L
    lq  <- log(1.3);   label("Intercompartmental clearance (Q, L/h)")                # Table 2: Q  = 1.3 L/h
    lvp <- log(10.2);  label("Peripheral volume of distribution (V2, L)")            # Table 2: V2 = 10.2 L

    # Inter-individual variability (Hodiamont 2017 Table 2 IIV block).
    # The Table 2 caption defines CV% = sqrt(exp(omega^2) - 1) * 100, so
    # the variance on the log scale is omega^2 = log(1 + CV^2).
    #   IIV CL : 75.0% CV  -> log(1 + 0.75^2)  = 0.446287
    #   IIV V1 : 27.0% CV  -> log(1 + 0.27^2)  = 0.070365
    # CL and V1 are correlated (r = 0.21 in Table 2); the lower-triangular
    # covariance is cov(CL, V1) = r * omega_CL * omega_V1
    #   = 0.21 * sqrt(0.446287) * sqrt(0.070365) = 0.037214.
    etalcl + etalvc ~ c(0.446287,
                        0.037214, 0.070365)   # Hodiamont 2017 Table 2: IIV CL 75.0%, IIV V1 27.0%, corr 0.21
    # Q and V2 had no IIV in the final model; the paper reports IIV terms
    # only for CL and V1 (Table 2 IIV block).

    # Residual error (Hodiamont 2017 Table 2 Residual variability block:
    # combined proportional + additive).
    propSd <- 0.194;  label("Proportional residual error (fraction)")   # Table 2: 19.4%
    addSd  <- 0.13;   label("Additive residual error (mg/L)")           # Table 2: 0.13 (mg/L; concentration units throughout the paper)
  })

  model({
    # Individual PK parameters. The final model retains no covariate
    # effects (Hodiamont 2017 Results "Model" paragraph); body weight and
    # other covariates were tested but did not improve fit. Q and V2 are
    # typical-value only because the paper does not report IIV on them.
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q  <- exp(lq)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment IV PK. Doses go directly to the central compartment
    # as a 30-min infusion in the source paper; the library model does
    # not hard-code the infusion duration so users can specify rate/dur
    # per dose in their event table.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Plasma concentration: dose in mg, vc in L -> mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
