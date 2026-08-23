Zhang_2025_fluconazole <- function() {
  description <- paste(
    "One-compartment population PK model for fluconazole in critically ill adults with acute renal",
    "failure receiving continuous renal replacement therapy (CRRT), with a tandem CRRT compartment",
    "from which extracorporeal clearance removes drug. Total clearance is the sum of residual body",
    "clearance (estimated separately in acute renal failure and in normal renal function) and the",
    "per-record CRRT clearance QEFF, gated on RRT_CRRT_ACTIVE so that clearance during an",
    "interruption reduces to body clearance alone. Body weight scales the central volume as a power",
    "function of WT/70. IV infusion only; no absorption. Fitted to 297 plasma concentrations",
    "digitised from six published CRRT case series.",
    sep = " "
  )
  reference <- paste(
    "Zhang S, Zhang W, Wu T, Qin Y, Pei Q. Optimizing fluconazole dosing in acute renal failure",
    "patients undergoing continuous renal replacement therapy: A population",
    "pharmacokinetic/pharmacodynamic study. Front Pharmacol. 2025;16:1564070.",
    "doi:10.3389/fphar.2025.1564070.",
    sep = " "
  )
  vignette <- "Zhang_2025_fluconazole"

  # The second state is the authors' "CRRT compartment" (Figure 2: "CLcrrt, CRRT compartment
  # clearance; Vcrrt, the apparent distribution volume of the CRRT compartment"). It is a
  # kinetic construct, not an anatomical space or a circuit volume: it is the lumped body space
  # accessible to the extracorporeal filter, from which CRRT removes drug at CLcrrt/Vcrrt.
  # Declared paper-specific rather than promoted to a canonical compartment per operator ruling
  # (sidecar request-001 question 1, answered 2026-08-20, option A) -- reversible if a second
  # CRRT paper needs the same role.
  paper_specific_compartments <- "crrt"

  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "fluconazole", units = "mg", specimen = "plasma", verified = TRUE),
    # specimen = "not applicable": the source defines this state as a kinetic compartment for
    # CRRT elimination (Results 3.2, "the CRRT compartment describes drug elimination via CRRT")
    # and never assigns it a biological matrix. Vcrrt = 23.5 L is a body-scale apparent volume,
    # so it is NOT the extracorporeal circuit (contrast Leuppi-Taegtmeyer_2019_colistin.R, whose
    # filter/cartridge states carry litre-scale priming volumes).
    crrt    = list(analyte = "fluconazole", units = "mg", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight, scaling the central compartment volume as (WT/70)^e_wt_vc.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Reference weight 70 kg, taken verbatim from the supplement's $PK block",
        "(V1 = THETA(3)*(WT/70)**THETA(6)*EXP(ETA(2))) and confirmed by the Discussion",
        "(\"the body weight-to-70 kg ratio was integrated into the Vc calculation using a power",
        "function\"). Cohort median 77 kg, range 48-272 kg (Table 1, Overall). Body weight was",
        "NOT retained on clearance.",
        sep = " "
      ),
      source_name        = "WT"
    ),
    DIS_ARF = list(
      description        = paste(
        "Binary acute-renal-failure indicator selecting which of the two residual body-clearance",
        "estimates applies: 1 = acute renal failure (lcl_arf, 0.41 L/h), 0 = normal renal function",
        "(lcl_nrf, 1.25 L/h).",
        sep = " "
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (normal renal function)",
      notes              = paste(
        "The supplement's $PK block selects the stratum on the urine-output category column:",
        "CL = THETA(1)*EXP(ETA(1)) by default, overridden by IF (URINE.EQ.4) CL = THETA(2)*EXP(ETA(1)).",
        "URINE == 4 is the normal-renal-function subject, so DIS_ARF = 1 - (URINE == 4). The paper",
        "reports the covariate as ARF Yes/No (Table 1: 15 ARF, 1 non-ARF -- a liver-transplant",
        "recipient) and names the two estimates CLbody_arf / CLbody_nrf (Table 2). Of the 15 renal-",
        "failure patients, seven were anuric and seven oliguric (Results 3.1). Note the single",
        "non-ARF subject means lcl_nrf rests on one patient's data.",
        sep = " "
      ),
      source_name        = "URINE"
    ),
    QEFF = list(
      description        = paste(
        "Per-record CRRT clearance of fluconazole (the paper's CLcrrt), supplied as a data column",
        "rather than estimated. Removes drug from the CRRT compartment at QEFF/vcrrt.",
        sep = " "
      ),
      units              = "L/h",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Computed per patient from CRRT mode and circuit flows (Supplementary Material Section 1):",
        "CVVH, CL_CRRT = Q_uf * S_c; CVVHD, CL_CRRT = Q_d * S_d; CVVHDF,",
        "CL_CRRT = (Q_uf + Q_d) * S_d, where Q_uf is the ultrafiltration flow rate (L/h), Q_d the",
        "dialysate flow rate (L/h), and S_c / S_d the sieving / saturation coefficient derived as",
        "AUC_filtrate / AUC_plasma. All three modes collapse to (Q_uf + Q_d) * S when the inactive",
        "flow is zero, which is how the paper's simulation tiers are specified: the Shiny user guide",
        "(Supplementary Material Section 4) defines \"the RRT dose\" as \"the sum of the dialysate flow",
        "rate (Q_D) and the ultrafiltration rate (Q_UF), mL/Kg/h\". With the cohort sieving/saturation",
        "coefficient 0.67 (Table 1, Overall), a CRRT dose in mL/kg/h converts as",
        "QEFF [L/h] = dose [mL/kg/h] * WT [kg] * 0.67 / 1000. Observed cohort range 8.0-64.8 mL/h/kg",
        "(Table 1). Carries no IIV: it is measured circuit hardware, not an estimated parameter.",
        sep = " "
      ),
      source_name        = "CLCRRT"
    ),
    RRT_CRRT_ACTIVE = list(
      description        = paste(
        "Time-varying gate for whether the continuous renal replacement therapy circuit is running.",
        "1 while running, 0 while interrupted or in a subject not on CRRT.",
        sep = " "
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CRRT not running)",
      notes              = paste(
        "Reproduces the supplement's $PK gate IF (CRRTYN.EQ.0) CL_crrt = 0 /",
        "IF (CRRTYN.EQ.1) CL_crrt = CLCRRT verbatim, so removal via the CRRT compartment ceases",
        "when the circuit stops and total clearance reduces to body clearance alone. Every subject",
        "in the model-building dataset was on CRRT (Table 1); the gate matters for the paper's",
        "simulations, which contrast on-CRRT against off-CRRT (the 0 mL/kg/h row of Tables 3 and 4)",
        "and for the intermittent-session scenarios the companion Shiny application supports.",
        sep = " "
      ),
      source_name        = "CRRTYN"
    )
  )

  # Screened during stepwise covariate modelling (Methods 2.2: "Potential covariates such as age,
  # body weight, and filter membrane area were considered continuous factors, while sex, ARF, and
  # filter membrane type were regarded as categorical variables") but NOT retained in the final
  # model, which kept only renal failure on clearance and body weight on volume (Results 3.2).
  # The corresponding $INPUT columns (AGE, GENDER, MEMAREA, MRMTYPE) are present in the dataset
  # but absent from the final $PK block. No point estimates are reported for any of them.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at the time of the modelled CRRT treatment.",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Screened as a continuous covariate and not retained. Cohort median 59.5 years, range",
        "32-82 years (Table 1, Overall); note the Results 3.1 text instead states a median of",
        "72 years, which disagrees with Table 1. Source column AGE.",
        sep = " "
      )
    ),
    SEXF = list(
      description = "Female sex indicator.",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened as a categorical covariate and not retained. Cohort 9 male / 7 female",
        "(Table 1, Overall). Source column GENDER; the paper does not state its coding, so the",
        "direction of the SEXF transformation could not be confirmed -- immaterial here because",
        "the covariate does not enter the final model.",
        sep = " "
      )
    ),
    MEMBRANE_AREA = list(
      description = "Surface area of the CRRT filter membrane.",
      units       = "m^2 (assumed; the paper does not state the unit)",
      type        = "continuous",
      notes       = paste(
        "Screened as a continuous covariate and not retained; no values are tabulated. Source",
        "column MEMAREA. NOT a canonical register name: no entry in",
        "inst/references/covariate-columns.md covers dialyser membrane surface area, and none was",
        "minted here because Zhang 2025 screened and rejected the covariate. A future model that",
        "RETAINS a dialyser-membrane covariate should file a naming sidecar rather than reuse this",
        "ad-hoc name.",
        sep = " "
      )
    ),
    MEMBRANE_TYPE = list(
      description = "CRRT filter membrane material / type.",
      units       = "(categorical)",
      type        = "categorical",
      notes       = paste(
        "Screened as a categorical covariate and not retained; the categories are not enumerated.",
        "Source column MRMTYPE. NOT a canonical register name, for the same reason given under",
        "MEMBRANE_AREA.",
        sep = " "
      )
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 16L,
    n_studies       = 6L,
    n_observations  = 297L,
    age_range       = "32-82 years",
    age_median      = "59.5 years",
    weight_range    = "48-272 kg",
    weight_median   = "77 kg",
    sex_female_pct  = 43.75,
    race_ethnicity  = NULL,
    disease_state   = paste(
      "Critically ill adults receiving continuous renal replacement therapy. 15 of 16 had acute",
      "renal failure (seven anuric, seven oliguric); the remaining subject was a liver-transplant",
      "recipient with normal renal function.",
      sep = " "
    ),
    dose_range      = "50-1200 mg fluconazole IV, infused over 0.5-6 h",
    regions         = "not reported",
    renal_function  = paste(
      "Acute renal failure in 15 of 16 subjects, all on CRRT. The single normal-renal-function",
      "subject is the sole source of information for lcl_nrf.",
      sep = " "
    ),
    crrt_modality   = paste(
      "CVVH / CVVHD / CVVHDF counts of 13 / 1 / 14 in the Table 1 Overall row. These sum to 28",
      "across 16 patients, so they enumerate treatment periods rather than patients (Valtonen 1997",
      "and Muhl 2000 each studied both CVVH and CVVHDF).",
      sep = " "
    ),
    crrt_parameters = paste(
      "Blood flow rate Qb 75-180 mL/min; CRRT dose 8.0-64.8 mL/h/kg; sieving / saturation",
      "coefficient Sc/Sd median 0.67 overall, per-study medians 0.49-0.87 (Table 1).",
      sep = " "
    ),
    notes           = paste(
      "Demographics are from Table 1. The dataset is not original patient-level data: 297 plasma",
      "concentration points were digitised with WebPlotDigitizer 4.3 from published",
      "concentration-time figures in six papers (Wolter 1994, Valtonen 1997, Muhl 2000,",
      "Kishino 2001, Yagasaki 2003, Lopez and Phillips 2014), with missing covariate values imputed",
      "by the mean or median. Methods 2.1 and Results 3.1 both say \"five\" publications, but Table 1",
      "tabulates six and the sixth (Kishino 2001) supplies the single normal-renal-function subject,",
      "so n_studies is recorded as 6. Digitisation accuracy was checked against a source with raw",
      "values available (Supplementary Table S1 / Figure S1; mean difference 0.0004 mg/L).",
      "Parameter uncertainty was obtained by sampling importance resampling (1000 proposal samples,",
      "1000 resamples) rather than bootstrap, because of the small number of subjects.",
      sep = " "
    )
  )

  ini({
    # All final estimates are from Table 2, column "Final PK model / Estimate". The supplement's
    # $THETA and $OMEGA blocks hold INITIAL values (0.4, 1.25, 42, 22, 36, 0.8 and 0.2, 0.1, 0.1,
    # 0.3) which differ from the final estimates and are deliberately not used here. Every
    # parameter carries an RSE% in Table 2, so none was fixed; nothing is wrapped in fixed().

    # Residual body clearance, estimated separately in each renal-function stratum within one
    # joint fit (supplement $PK: CL = THETA(1)*EXP(ETA(1)), overridden by
    # IF (URINE.EQ.4) CL = THETA(2)*EXP(ETA(1))). Both strata share the single IIV term etalcl,
    # so per the stratum-suffix convention only the typical values are suffixed.
    lcl_arf <- log(0.41);   label("Residual body clearance, acute renal failure (L/h)")      # Table 2, CLbody_arf = 0.41 (RSE 15%, SIR 95% CI 0.30-0.52)
    lcl_nrf <- log(1.25);   label("Residual body clearance, normal renal function (L/h)")    # Table 2, CLbody_nrf = 1.25 (RSE 23%, SIR 95% CI 0.69-1.81)

    lvc     <- log(37.90);  label("Central compartment volume at 70 kg (L)")                 # Table 2, Vc = 37.90 (RSE 11%, SIR 95% CI 30.65-45.27)
    lvcrrt  <- log(23.50);  label("CRRT compartment volume (L)")                             # Table 2, Vcrrt = 23.50 (RSE 12%, SIR 95% CI 18.21-28.82)
    lq      <- log(38.80);  label("Intercompartmental clearance between central and CRRT compartments (L/h)")  # Table 2, Q = 38.80 (RSE 25%, SIR 95% CI 20.07-48.95)

    e_wt_vc <- 0.799;       label("Power exponent on (WT/70) for central volume (unitless)")  # Table 2, "WT on Vc" = 0.799 (RSE 7%, SIR 95% CI 0.690-0.914); supplement $PK THETA(6)

    # IIV. Table 2 heads this block "Interindividual variability (omega^2)", so the tabulated
    # values are variances on the log scale and are used directly. Approximate CVs from
    # sqrt(exp(omega^2) - 1): 45%, 32%, 30%, 71%.
    etalcl    ~ 0.187;      label("IIV variance on log residual body clearance (both strata share ETA(1))")  # Table 2, IIV CL = 0.187 (RSE 18%, shrinkage 15.9%)
    etalvc    ~ 0.097;      label("IIV variance on log central volume")                      # Table 2, IIV Vc = 0.097 (RSE 29%, shrinkage 7.8%)
    etalvcrrt ~ 0.084;      label("IIV variance on log CRRT compartment volume")             # Table 2, IIV Vcrrt = 0.084 (RSE 59%, shrinkage 36.4%)
    etalq     ~ 0.412;      label("IIV variance on log intercompartmental clearance")        # Table 2, IIV Q = 0.412 (RSE 35%, shrinkage 44.5%)

    # Residual error. The supplement's $ERROR block is Y = F + EPS(1) -- pure additive, with no
    # proportional term. Table 2 heads the block "Residual variability (sigma)" and annotates the
    # row "Additive error (mg/L)"; a variance would carry mg^2/L^2, so 0.375 is read as an SD.
    addSd <- 0.375;         label("Additive residual error (mg/L)")                          # Table 2, Additive error = 0.375 (RSE 13%, SIR 95% CI 0.298-0.473)
  })

  model({
    # Residual body clearance: the two stratum estimates selected by DIS_ARF, sharing one IIV
    # term. Written as a single expression so that DIS_ARF = 1 gives exp(lcl_arf + etalcl) and
    # DIS_ARF = 0 gives exp(lcl_nrf + etalcl), reproducing the supplement's $PK IF-block exactly.
    cl     <- exp(lcl_arf * DIS_ARF + lcl_nrf * (1 - DIS_ARF) + etalcl)

    vc     <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    vcrrt  <- exp(lvcrrt + etalvcrrt)
    q      <- exp(lq + etalq)

    # Extracorporeal clearance, supplied per record and switched off when the circuit stops
    # (supplement $PK: IF (CRRTYN.EQ.0) CL_crrt = 0; IF (CRRTYN.EQ.1) CL_crrt = CLCRRT).
    clcrrt <- QEFF * RRT_CRRT_ACTIVE

    # Micro-constants, named as in the supplement's $PK block (K10, K12, K21, K20). Note that
    # K20 removes drug from the CRRT compartment, NOT from central -- that extra elimination
    # path is what distinguishes this structure from a conventional two-compartment model.
    kel    <- cl / vc            # K10, body elimination from central (1/h)
    k12    <- q / vc             # K12, central -> CRRT compartment (1/h)
    k21    <- q / vcrrt          # K21, CRRT compartment -> central (1/h)
    k20    <- clcrrt / vcrrt     # K20, extracorporeal removal from the CRRT compartment (1/h)

    # Supplement $DES, verbatim:
    #   DADT(1) = -K10*A(1) - K12*A(1) + K21*A(2)
    #   DADT(2) =  K12*A(1) - K20*A(2) - K21*A(2)
    d/dt(central) <- -kel * central - k12 * central + k21 * crrt
    d/dt(crrt)    <-  k12 * central - k20 * crrt    - k21 * crrt

    # Supplement $MODEL declares COMP=(CENTRAL, DEFDOS, DEFOBS) with S1 = V1, so plasma
    # concentration is the central amount over the central volume, in mg/L.
    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
