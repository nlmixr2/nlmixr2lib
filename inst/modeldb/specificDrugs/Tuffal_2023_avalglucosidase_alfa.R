Tuffal_2023_avalglucosidase_alfa <- function() {
  description <- paste(
    "Concatenated 3-compartment population PK model for avalglucosidase",
    "alfa (enzyme replacement therapy) in adolescent and adult patients",
    "with late-onset Pompe disease (Tuffal 2023). Pooled analysis of 2042",
    "plasma drug-activity determinations from 75 patients across three",
    "trials (phase I/II NCT01898364, its follow-up NCT02032524, and phase",
    "III NCT02782741) at 5, 10, and 20 mg/kg IV Q2W. Elimination from the",
    "central compartment is the sum of a linear clearance (CL) and a",
    "parallel saturable Michaelis-Menten clearance (Vmax, Km). The two",
    "peripheral compartments are arranged in SERIES rather than in",
    "parallel: central exchanges bidirectionally with peripheral1 (Q2),",
    "peripheral1 feeds peripheral2 one way (Q3), and peripheral2 returns",
    "drug directly to central one way (Qpc). That cycle is what produces",
    "the two chronological kinetic sequences seen in the data - a first",
    "exposure-driving sequence contributing about 99% of AUC, then a slow",
    "low-concentration rebound peaking near 87 h and remaining above the",
    "LLOQ out to the next dose at 336 h. The authors explicitly tested and",
    "rejected a bidirectional Q3 ('Markedly increased OFV'). Q2, V2, Q3,",
    "and V3 were fixed after preliminary screening to avoid identifiability",
    "problems. No covariates are retained: 12 demographic / laboratory",
    "covariates and 3 antidrug-antibody parameterisations were screened and",
    "none qualified, so the screened set is documented in",
    "covariatesDataExcluded rather than covariateData. Allometric scaling",
    "was tested and not retained, so the model carries no body-weight term."
  )
  reference <- paste(
    "Tuffal G, Tiraboschi G, Hurbin F, Boittet P, Palmer R,",
    "Martinez J-M, Fabre D. Population Pharmacokinetic Modeling and",
    "Determination of Individual Exposure to Avalglucosidase Alfa in",
    "Adolescent and Adult Patients With Late-Onset Pompe Disease:",
    "Analysis of Pooled Data From Phase I to III Clinical Trials.",
    "Ther Drug Monit. 2023;45(5):644-652.",
    "doi:10.1097/FTD.0000000000001086."
  )
  vignette <- "Tuffal_2023_avalglucosidase_alfa"

  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "ug/mL"
  )

  # No covariates are retained in the final model. Every covariate the paper
  # screened is documented in covariatesDataExcluded below.
  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Derived mechanically; verified = FALSE means it has
  # NOT been checked against the source paper.
  compartmentData <- list(
    central     = list(analyte = "avalglucosidase alfa", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "avalglucosidase alfa", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral2 = list(analyte = "avalglucosidase alfa", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list()

  covariatesDataExcluded <- list(
    WT = list(
      description        = "Baseline body weight. Screened as a continuous covariate on CL and V1, and separately as the scaling variable for allometric scaling of central clearance and volume (Methods, 'Population Pharmacokinetic Modeling' item 1).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Survived backward deletion on CL IIV (removal raised OFV by more than the 10.84 criterion) but was NOT retained: the resulting model failed the covariance-step acceptance criteria (correlated parameters and/or several RSE much greater than 50%), the effect on CL IIV was not significant, and the single-covariate bootstrap qualification failed (Results, 'Population Pharmacokinetic Model'). Allometric scaling was likewise tested and not retained, so the final model contains no body-weight term at all. Reported exposure nevertheless tracks weight: median AUC0-336h was 32% lower below 50 kg and 41% higher at or above 100 kg relative to 50-100 kg (Results, 'Patient Exposure'; Supplemental Digital Content 6). Population mean weight 75.9 kg (SD 20.1), Table 1."
    ),
    AGE = list(
      description        = "Baseline age. Screened as a continuous demographic covariate (Methods).",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Not retained. Population mean 46.0 years (SD 15.1), Table 1; a single patient was under 18 years."
    ),
    SEXF = list(
      description        = "Sex indicator. Screened as a categorical demographic covariate (Methods).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Not retained. Paper reports counts rather than an indicator column: 36 of 75 female (48%), Table 1."
    ),
    RACE_WHITE = list(
      description        = "White / Caucasian race indicator. Part of the 'ethnicity' categorical covariate screened in the Methods.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Caucasian)",
      notes              = "Not retained. Table 1 reports Caucasian 68 of 75 (90.7%). The paper screened ethnicity as a single categorical covariate; it is decomposed here into the canonical one-hot indicators so the screened groups are individually documented. Table 3 labels the same variable 'Phenotype'."
    ),
    RACE_BLACK = list(
      description        = "Black race indicator. Part of the 'ethnicity' categorical covariate screened in the Methods.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not Black)",
      notes              = "Not retained. Table 1 reports Black 2 of 75 (2.67%)."
    ),
    RACE_ASIAN = list(
      description        = "Asian race indicator. Part of the 'ethnicity' categorical covariate screened in the Methods.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not Asian)",
      notes              = "Not retained. Table 1 reports Asian 3 of 75 (4.00%); a further 2 patients (2.67%) were recorded as 'Other'."
    ),
    CRCL = list(
      description        = "Baseline renal function. The paper screened renal function under two parameterisations: creatinine clearance, and glomerular filtration rate estimated by the Modification of Diet in Renal Disease (MDRD) formula (Methods).",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Neither parameterisation was retained. Both map to the single canonical renal-function column, which covers creatinine-based estimates and measured GFR alike; the two parameterisations are recorded here rather than split across two columns because neither entered the final model. Table 3 stratifies exposure by GFR: 64 of 70 patients at or above 90 mL/min, 6 in the 60-89 range."
    ),
    CPK = list(
      description        = "Baseline creatine kinase. Screened as a continuous laboratory covariate (Methods).",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Not retained. Baseline value used, per the Methods convention 'Unless otherwise stated, the baseline value was used as the covariate.' Paper calls this 'creatine kinase'; units are not stated."
    ),
    ALP = list(
      description        = "Baseline alkaline phosphatase. Screened as a continuous laboratory covariate (Methods).",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Not retained. Units not stated in the paper."
    ),
    ALT = list(
      description        = "Baseline alanine aminotransferase. Screened as a continuous laboratory covariate (Methods).",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Not retained. Units not stated in the paper."
    ),
    AST = list(
      description        = "Baseline aspartate aminotransferase. Screened as a continuous laboratory covariate (Methods).",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Not retained, despite surviving backward deletion on CL IIV alongside body weight (removal raised OFV by more than the 10.84 criterion). Rejected for the same reasons as WT: covariance-step acceptance criteria not met, effect on CL IIV not significant, and single-covariate bootstrap qualification failed (Results, 'Population Pharmacokinetic Model'). Units not stated in the paper."
    ),
    TBILI = list(
      description        = "Baseline total bilirubin. Screened as a continuous laboratory covariate (Methods).",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Not retained. The paper says only 'bilirubin' and states no units; total bilirubin is assumed."
    ),
    ALB = list(
      description        = "Baseline albumin. Screened as a continuous laboratory covariate (Methods).",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Not retained. Units not stated in the Methods, but Table 3 stratifies exposure at a 45 g/L cut point (28 patients below, 42 at or above), which implies g/L."
    ),
    ADA_POS = list(
      description        = "Antidrug-antibody positivity. The paper screened two of its three ADA parameterisations through this concept: a subject-level categorical covariate (ADA occurring at least once versus never during follow-up) and a longitudinal categorical covariate (positive or not at the current time).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ADA-negative)",
      notes              = "Neither parameterisation was retained. Both are recorded against this single canonical column: the ever-positive reading is time-fixed and the longitudinal reading is time-varying, but the column encoding is identical and neither entered the final model. See Methods, 'Population Pharmacokinetic Modeling', and Supplemental Digital Content 2 for the ADA assays."
    ),
    ADA_TITER = list(
      description        = "Continuous antidrug-antibody titer. The third of the paper's three ADA parameterisations (Methods).",
      units              = "titer (reciprocal dilution; assay described in Supplemental Digital Content 2)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Not retained. Qualitative and semiquantitative titer assays were performed for the phase III study only; two further assays characterised the neutralising ability of ADA (enzymatic-activity inhibition and uptake inhibition). Table 3 does not stratify exposure by ADA status."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 75L,
    n_studies      = 3L,
    age_range      = "at least 3 years by protocol eligibility; the analysis population was adolescent and adult, with 1 patient under 18 years, 59 aged 18-64, and 10 aged 65 or older among the 70 who received 20 mg/kg (Table 3)",
    age_median     = "mean 46.0 years (SD 15.1) overall; 46.0 (16.6) in phase I/II and 46.1 (14.5) in phase III (Table 1)",
    weight_range   = "not reported as a range; exposure was stratified below 50 kg (5 patients), 50-99 kg (55), and at or above 100 kg (10) among the 70 who received 20 mg/kg (Table 3)",
    weight_median  = "mean 75.9 kg (SD 20.1) overall; 72.0 (14.5) in phase I/II and 77.8 (22.1) in phase III (Table 1)",
    sex_female_pct = 48,
    race_ethnicity = c(Caucasian = 90.7, Black = 2.67, Asian = 4.00, Other = 2.67),
    disease_state  = "late-onset Pompe disease (LOPD), confirmed by GAA enzyme deficiency from any tissue source, 2 confirmed pathogenic GAA variants, or both; 14 of 75 (18.7%) had been pretreated with alglucosidase alfa for at least 9 months and 61 (81.3%) were treatment-naive",
    dose_range     = "5, 10, or 20 mg/kg IV Q2W. Each infusion was given stepwise to limit hypersensitivity reactions: 1 mg/kg/h for 30 min, 3 mg/kg/h for 30 min, 5 mg/kg/h for 30 min, then 7 mg/kg/h until the planned dose was delivered. 70 of the 75 patients received 20 mg/kg (51 from phase III plus 19 from phase I/II who started at or shifted to the top dose).",
    regions        = "not reported; multicentre across NCT01898364, NCT02032524, and NCT02782741",
    notes          = "Baseline demographics are in Table 1; the per-study design, dosing, and sampling schedules are in Supplemental Digital Content 1 Table 1. Data attrition: 2374 concentrations collected, 2066 above the LLOQ, 9 removed before modelling (3 predoses above LLOQ, 1 aberrant trough, 5 from a patient with truncated prior-dose information) giving 2057, then 15 conditional-weighted-residual outliers (0.7%) removed, giving the final 2042 used for fitting. Concentrations are plasma enzymatic activity, assayed fluorometrically over a validated range of 0.0125 (LLOQ) to 3.0 ug/mL. Fit in NONMEM 7.4.1 with FOCE-I. Final-model performance: mean bias -2.66% and RMSE 30.7% for IPRED, -0.433% and 38.9% for PRED; AAFE 1.36 (IPRED) and 1.72 (PRED). Bootstrap convergence 86% of 1000 runs."
  )

  ini({
    # Structural parameters -- Tuffal 2023 Table 2, "Fixed effects" block.
    # Estimated parameters (RSE reported):
    lcl   <- log(0.869);  label("Linear clearance from the central compartment CL (L/h)")            # Table 2: CL = 0.869 L/h (RSE 0.09%, 95% CI 0.868-0.871; bootstrap median 0.884)
    lvc   <- log(3.35);   label("Central compartment volume V1 (L)")                                 # Table 2: V1 = 3.35 L (RSE 0.07%, 95% CI 3.35-3.36; bootstrap median 3.44)
    lvmax <- log(9.82);   label("Michaelis-Menten maximum elimination rate Vmax (mg/h)")             # Table 2: Vmax = 9.82 mg/h (RSE 0.19%, 95% CI 9.79-9.86; bootstrap median 10.6)
    lkm   <- log(0.451);  label("Michaelis-Menten half-saturation constant Km (ug/mL)")              # Table 2: Km = 0.451 ug/mL (RSE 0.12%, 95% CI 0.450-0.452; bootstrap median 0.493)
    lqpc  <- log(0.0206); label("Back-redistribution clearance peripheral2 -> central Qpc (L/h)")    # Table 2: Qpc = 0.0206 L/h (RSE 0.330%, 95% CI 0.0204-0.0207; bootstrap median 0.0183)

    # Fixed parameters. Results, "Population Pharmacokinetic Model": "To avoid
    # identifiability issues in the model, Q2, Q3, V2, and V3 values were fixed
    # after the preliminary model screening and sensitivity analysis." Table 2
    # marks each of these four rows "Fixed value" with no RSE, CI, or bootstrap.
    lq    <- fixed(log(0.254)); label("Inter-compartmental clearance central <-> peripheral1 Q2 (L/h)")  # Table 2: Q2 = 0.254 L/h (Fixed value)
    lvp   <- fixed(log(296));   label("First peripheral compartment volume V2 (L)")                      # Table 2: V2 = 296 L (Fixed value)
    lq2   <- fixed(log(1.87));  label("One-way inter-compartmental clearance peripheral1 -> peripheral2 Q3 (L/h)")  # Table 2: Q3 = 1.87 L/h (Fixed value)
    lvp2  <- fixed(log(1.31));  label("Second peripheral compartment volume V3 (L)")                     # Table 2: V3 = 1.31 L (Fixed value)

    # IIV -- Tuffal 2023 Table 2, "Inter-individual variability" block. The
    # Estimate column carries the log-scale variance omega^2 and the
    # parenthesised value is the corresponding CV%. Confirmed against the
    # paper's own printed CV% for all five etas via CV = sqrt(exp(omega^2) - 1):
    #   0.0870 -> 30.2%   0.0160 -> 12.7%   0.125 -> 36.6%
    #   0.171  -> 43.1%   1.87   -> 234%
    # matching Table 2 exactly, and matching the Results sentence "IIV in model
    # parameters ranged from 12.7% to 43.1% for CL, V1, Vmax, and Km and up to
    # 234% for Qpc". Results also states "A diagonal matrix was used to explain
    # IIV", so there are no off-diagonal terms. "The log-normal IIV of the
    # parameters was used."
    etalcl   ~ 0.0870   # Table 2: omega^2 CL   = 0.0870 (CV 30.2%), RSE 20.5%
    etalvc   ~ 0.0160   # Table 2: omega^2 V1   = 0.0160 (CV 12.7%), RSE 32.1%
    etalvmax ~ 0.125    # Table 2: omega^2 Vmax = 0.125  (CV 36.6%), RSE 34.8%
    etalkm   ~ 0.171    # Table 2: omega^2 Km   = 0.171  (CV 43.1%), RSE 10.1%
    etalqpc  ~ 1.87     # Table 2: omega^2 Qpc  = 1.87   (CV 234%),  RSE 1.21%

    # Residual error. Results: "the residual variability was best described
    # using a proportional error model". Table 2 reports sigma^2 = 0.125
    # (CV 35.3%); for a proportional model the SD is sqrt(0.125) = 0.3536,
    # which reproduces the printed 35.3% exactly. Corroborated by the
    # reported RMSE of 30.7% (IPRED) and 38.9% (PRED), which bracket it.
    propSd <- 0.354; label("Proportional residual error (fraction)")  # Table 2: residual variability sigma^2 = 0.125 (CV 35.3%), RSE 0.01% -> propSd = sqrt(0.125)
  })

  model({
    # Individual parameters. Estimated parameters carry log-normal IIV; the
    # four parameters Table 2 marks "Fixed value" carry none.
    cl   <- exp(lcl + etalcl)
    vc   <- exp(lvc + etalvc)
    vmax <- exp(lvmax + etalvmax)
    km   <- exp(lkm + etalkm)
    qpc  <- exp(lqpc + etalqpc)
    q    <- exp(lq)
    vp   <- exp(lvp)
    q2   <- exp(lq2)
    vp2  <- exp(lvp2)

    # Concentrations driving the clearance terms. Amounts are in mg and volumes
    # in L, so every concentration is in mg/L = ug/mL -- the units of Km
    # (0.451 ug/mL) and of the reported observations.
    Cc  <- central / vc
    Cp1 <- peripheral1 / vp
    Cp2 <- peripheral2 / vp2

    # Concatenated ("series") 3-compartment system with back-redistribution,
    # per Tuffal 2023 Figure 1 inset and the Results parameter list "CL, Vmax,
    # and Km, describing the 2 parallel clearances from the central
    # compartment, and Q2, Q3, Qpc, V1, V2, and V3, describing drug
    # distribution."
    #
    #   central <--Q2--> peripheral1 --Q3--> peripheral2 --Qpc--> central
    #      |  linear CL (out)
    #      |  Michaelis-Menten Vmax/Km (out)
    #
    # Only Q2 is bidirectional. Q3 and Qpc are one-way, which is what closes
    # the loop and generates the second kinetic sequence. The one-way reading is
    # the authors' own: Supplemental Digital Content 1 Table 2 compares a
    # "Q3 one way" run against a "Q3 two ways" run and rejects the latter with
    # the comment "Markedly increased OFV"; every later development step,
    # including "Final model (No allometry)", descends from the one-way branch.
    # Elimination is from the central compartment only, and both clearances act
    # in parallel there.
    d/dt(central) <- -cl * Cc -
                      vmax * Cc / (km + Cc) -
                      q * Cc + q * Cp1 +
                      qpc * Cp2
    d/dt(peripheral1) <- q * Cc - q * Cp1 - q2 * Cp1
    d/dt(peripheral2) <- q2 * Cp1 - qpc * Cp2

    Cc ~ prop(propSd)
  })
}
