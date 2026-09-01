ButraguenoLaiseca_2022_piperacillin <- function() {
  description <- paste(
    "Two-compartment population PK model for intravenous piperacillin in 32",
    "critically ill children in a paediatric ICU, 13 of whom were receiving",
    "continuous kidney replacement therapy (CKRT) in continuous venovenous",
    "haemodiafiltration modality (Butragueno-Laiseca 2022). Plasma,",
    "pre-filter, post-filter, effluent and urine concentrations were fitted",
    "simultaneously, which is what identifies the three parallel elimination",
    "pathways. Elimination from the central compartment is the ADDITIVE sum",
    "of a renal arm (CLR, 1.3 L/h typical, encoded as lcl_renal, scaled",
    "linearly by estimated glomerular filtration rate and by height), a",
    "non-renal arm (CLM, 0.5 L/h typical, encoded as lcl_nonren, scaled",
    "linearly by body weight) and a haemofilter arm (CLKRT, 1.34 L/h typical",
    "for the large filter, encoded as lcl_hemodialysis). The renal arm is",
    "gated off in patients without diuresis and loses its eGFR and height",
    "covariates in CKRT patients; the haemofilter arm is gated off when CKRT",
    "is not running. Haemofilter surface area is a categorical covariate on",
    "the CKRT arm with the LARGE 1.2 m2 filter as the reference (medium",
    "0.6 m2 multiplier 0.74; small 0.2 m2 multiplier 0.28). Both volumes",
    "scale linearly with body weight. Besides the plasma concentration Cc,",
    "the model returns three flow-normalised observables, each with its own",
    "residual error: the post-filter concentration (pre-filter concentration",
    "reduced by the single-pass filter extraction ratio), the effluent",
    "concentration and the urine concentration. Age, serum albumin,",
    "haematocrit and the other laboratory values were screened and not",
    "retained."
  )
  reference <- paste(
    "Butragueno-Laiseca L, Marco-Arino N, Troconiz IF, Grau S, Campillo N,",
    "Garcia X, Padilla B, Fernandez SN, Slocker M, Santiago MJ.",
    "Population pharmacokinetics of piperacillin in critically ill children",
    "including those undergoing continuous kidney replacement therapy.",
    "Clin Microbiol Infect. 2022;28(9):1287.e9-1287.e15.",
    "doi:10.1016/j.cmi.2022.03.031"
  )
  vignette <- "ButraguenoLaiseca_2022_piperacillin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central     = list(analyte = "piperacillin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "piperacillin", units = "mg", specimen = "plasma", verified = TRUE),
    urine       = list(analyte = "piperacillin", units = "mg", specimen = "urine",  verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters the non-renal clearance arm and both volumes of distribution as a plain",
        "LINEAR ratio (WGT/8.1), with NO allometric exponent -- Butragueno-Laiseca 2022",
        "Table 2 prints the Parameter model column as theta_CLM * WGT/8.1 and theta_V *",
        "WGT/8.1, and the Results are explicit that an exponent was tested and rejected:",
        "'Adding the allometric scaling factor to the CLM vs Body weight relationship did",
        "not improve the fit (p > 0.05), and fixing it to 0.75 increased the -2LL value.'",
        "The reference weight of 8.1 kg is stated in Table 2 footnote d ('Median body",
        "weight = 8.1 kg') and equals the cohort median (Results, Demographics: median",
        "(range) weight 8.1 kg (4 to 63 kg)). Weight was NOT retained on the renal",
        "clearance arm or on the inter-compartmental clearance (Results: body weight",
        "correlated with V1, V2 and CLM 'but not with CLR or inter-compartmental",
        "distribution clearance (CLD; p > 0.05)'). Treated as a time-fixed baseline value",
        "in the source analysis. Weight also acts as an implicit surrogate for haemofilter",
        "size, because filter size was assigned by weight band in the target-attainment",
        "simulations (0.2 m2 for 3-10 kg, 0.6 m2 for 10-30 kg, 1.2 m2 for 30-60 kg;",
        "Methods, Probability of target attainment)."
      ),
      source_name        = "WGT"
    ),
    HT = list(
      description        = "Body height",
      units              = "cm",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Multiplies the renal clearance arm as the plain linear ratio (HGT/69), with no",
        "exponent (Butragueno-Laiseca 2022 Table 2: CLR = theta_CLR * (eGFR/119.3) *",
        "(HGT/69)). The reference height of 69 cm is stated in Table 2 footnote c ('Median",
        "height = 69 cm') and is restated in the Results ('eGFR, height, and weight equal",
        "to the median of the studied population (8.1 kg, 69 cm, and 119.3 mL/min/1.73m2,",
        "respectively)'). The authors read height as a maturation marker rather than a size",
        "descriptor: 'our data supported the inclusion of height, which can be interpreted",
        "as a marker of organ maturation' (Discussion), offered as the reason no",
        "postmenstrual-age effect on clearance was found in a cohort whose youngest patient",
        "was three months old. Applies only to patients not undergoing CKRT -- see",
        "RRT_CRRT_ACTIVE. Table 1 reports height as mean (SD) 79.5 (30.7) cm without CKRT",
        "and 101 (33.2) cm with CKRT, and the two groups did not differ significantly."
      ),
      source_name        = "HGT"
    ),
    CRCL = list(
      description        = "Estimated glomerular filtration rate (Schwartz equation), BSA-normalised",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Multiplies the renal clearance arm as the plain linear ratio (eGFR/119.3), with no",
        "exponent (Butragueno-Laiseca 2022 Table 2). Estimated with the Schwartz formula",
        "(Methods, Pharmacokinetic analysis; the cited reference is Schwartz 1976,",
        "Pediatrics 58:259-263). The reference value of 119.3 mL/min/1.73 m2 is stated in",
        "Table 2 footnote b and is the cohort median for the non-CKRT group (Results,",
        "Demographics). The LINEAR (exponent-1) reading is confirmed arithmetically by the",
        "Discussion: 'an increase in eGFR of 10 mL/min/1.73m2 over the median value is",
        "associated with an increase in CLR and CL of 8% and 5.5%, respectively' --",
        "10/119.3 = 8.4%, which reproduces the 8% only if the ratio enters to the first",
        "power, and 1.3 * 0.084 / 1.8 = 6.1% reproduces the clearance figure. eGFR was",
        "measured only in patients WITHOUT CKRT (Table 1 shows a dash for the CKRT group)",
        "and the covariate is correspondingly gated off for CKRT patients -- Methods: 'The",
        "impact of the GFR ... was explored in the renal clearance (CLR) for those patients",
        "without CKRT', and Table 2 footnote a: 'In CKRT patients showing residual diuresis,",
        "CLR is given solely by the typical estimate theta_CLR.' Table 1 reports mean (SD)",
        "eGFR of 123 (55) mL/min/1.73 m2 in the non-CKRT group."
      ),
      source_name        = "eGFR"
    ),
    RRT_CRRT_ACTIVE = list(
      description        = "CKRT-active indicator (1 while continuous venovenous haemodiafiltration is running, 0 otherwise)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no CKRT running)",
      notes              = paste(
        "Does two separate jobs in this model. (1) It gates the haemofilter clearance arm",
        "ON: the supplementary material states that 'CLRenal and CLRRT were absent in",
        "patients without diuresis or without hemofilter'. (2) It gates the eGFR and height",
        "covariates on the RENAL arm OFF, because eGFR was neither measured nor modelled in",
        "CKRT patients -- Table 2 footnote a: 'In CKRT patients showing residual diuresis,",
        "CLR is given solely by the typical estimate theta_CLR.' Unlike the sibling",
        "ButraguenoLaiseca_2025_teicoplanin model, the renal arm is NOT set to zero on CKRT",
        "here: piperacillin CKRT patients with residual diuresis retain the typical renal",
        "clearance, and only patients without diuresis lose the arm entirely (that second,",
        "independent gate is carried by URINE_FLOW). The modality is continuous venovenous",
        "haemodiafiltration (Results, Demographics), hence CRRT rather than the",
        "intermittent-haemodialysis member of the RRT_<modality>_<kind> family. The ACTIVE",
        "rather than the STATUS member is used for consistency with the sibling model and",
        "because the quantity is physically time-varying (supplementary Table 3 reports",
        "filter running times of 10.5 to 60.5 h, so a filter change interrupts the circuit",
        "within a 10-day simulated course), although the source analysis treats CKRT",
        "membership as fixed per patient."
      ),
      source_name        = "Patient type (with or without CKRT)"
    ),
    FILT_SA_MED = list(
      description        = "Medium haemofilter indicator (0.6 m2 membrane surface area)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (the large 1.2 m2 reference filter, the small 0.2 m2 filter, or no CKRT)",
      notes              = paste(
        "One of the three levels of the haemofilter-surface-area covariate on the CKRT",
        "clearance arm (Butragueno-Laiseca 2022 Table 2: CLCKRT = theta_CLCKRT * theta_FILT,",
        "with theta_FILT_High = 1 Fixed, theta_FILT_Med = 0.74 (0.6-0.9) and theta_FILT_Low",
        "= 0.28 (0.2-0.4)). NOTE THE REFERENCE LEVEL: this paper fixes the LARGE 1.2 m2",
        "filter to 1 and estimates the medium and small multipliers BELOW it, which is the",
        "opposite direction from the sibling ButraguenoLaiseca_2025_teicoplanin model, whose",
        "reference is the small filter. The canonical column semantics are unchanged --",
        "FILT_SA_MED and FILT_SA_LARGE both 0 still denotes the SMALL filter -- and the",
        "paper's choice of reference level is absorbed into the effect expression in",
        "model(), where the small-filter indicator is derived as",
        "(1 - FILT_SA_MED - FILT_SA_LARGE). Reproduces the Results text exactly: 'The",
        "estimated CLCKRT for a filter surface of 1.2 m2 was 1.34 L/h ... For filter",
        "surfaces of 0.6 and 0.2 m2, CLCKRT showed estimates of 1.01 and 0.38 L/h' (1.34 *",
        "0.74 = 0.99 and 1.34 * 0.28 = 0.38, the first differing from the printed 1.01 only",
        "by the rounding of the two-decimal multiplier). Meaningful only when",
        "RRT_CRRT_ACTIVE = 1. Mutually exclusive with FILT_SA_LARGE. Supplementary Table 3",
        "records 6 small, 4 medium and 3 large filters among the 13 CKRT patients."
      ),
      source_name        = "FILT_Med"
    ),
    FILT_SA_LARGE = list(
      description        = "Large haemofilter indicator (1.2 m2 membrane surface area)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (the small 0.2 m2 filter, the medium 0.6 m2 filter, or no CKRT)",
      notes              = paste(
        "The REFERENCE level of the haemofilter-surface-area covariate in this paper",
        "(Butragueno-Laiseca 2022 Table 2: theta_FILT_High = 1 Fixed, where FILT_high is",
        "defined in the table footnote as 'filter surface 1.2 m2 (reference)'), so it",
        "carries no estimated multiplier of its own in ini(). The indicator is nonetheless",
        "load-bearing in model(): it is what distinguishes the reference large filter, whose",
        "multiplier is 1, from the small filter, which is what FILT_SA_MED and FILT_SA_LARGE",
        "both being 0 denotes and which carries the 0.28 multiplier. Meaningful only when",
        "RRT_CRRT_ACTIVE = 1. Mutually exclusive with FILT_SA_MED. Supplementary Table 3",
        "records 3 large filters among the 13 CKRT patients, with the shortest mean running",
        "time of the three sizes (10.5 h, SD 13.4)."
      ),
      source_name        = "FILT_high"
    ),
    BFR = list(
      description        = "Blood flow rate through the CKRT extracorporeal circuit",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters ONLY the post-filter observation equation, not the clearance model. The",
        "supplementary material gives CPost = CPre * (1 - CLRRT / phi_Pl,corr) (eq. 4),",
        "'where phi_Pl,corr and phi_Effl are the corrected plasma and total effluent flows.",
        "The value of phi_Pl,corr was calculated as phi_Blood x BPR, being phi_Blood the",
        "blood flow rate and BPR the blood to plasma concentration ratio, an additional",
        "parameter to be estimated from the model.' The main-text Results then pin BPR: 'The",
        "parameter blood to plasma ratio was not significantly different from 1 (p < 0.01)',",
        "so phi_Pl,corr = phi_Blood and -- unlike the sibling teicoplanin model, which",
        "encodes the correction as a plasma fraction 1 - HCT/100 -- NO haematocrit covariate",
        "is needed here. UNITS: supplementary Table 3 reports blood flow as an absolute",
        "70 (SD 31) mL/min. The main-text Results print the same quantity weight-normalised",
        "as '5 (SD 2) ml/kg/h', which is impossible: at 5 mL/kg/h a median 8.1 kg patient",
        "has a circuit flow of 0.041 L/h, less than a ninth of the 0.38 L/h CLCKRT of the",
        "filter such a patient would receive, so eq. 4 would return a negative post-filter",
        "concentration. The per-MINUTE reading reproduces the paper's own published",
        "extraction ratio; see the model file's post-filter block and the vignette Errata.",
        "Constant per subject during the study. Meaningful only when RRT_CRRT_ACTIVE = 1;",
        "converted to L/h inside model()."
      ),
      source_name        = "Blood flow"
    ),
    RRT_CRRT_EFFLUENT_FLOW = list(
      description        = "Total effluent flow leaving the haemofilter (dialysate + substitution fluid + net ultrafiltration)",
      units              = "mL/h",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters ONLY the effluent observation equation. Supplementary material eq. 5:",
        "CEffl = CLRRT * CPre / phi_Effl, with phi_Effl the 'total effluent flow' and 'The",
        "values of phi_Blood and phi_Effluent were measured during the course of the study.'",
        "Supplementary Table 3 reports it weight-normalised by filter size -- 78 (SD 24)",
        "mL/kg/h for the 0.2 m2 filter, 54 (5) for 0.6 m2 and 60 (15) for 1.2 m2 -- and the",
        "main-text Results give the cohort median as 68 (SD 20) mL/kg/h; multiply by body",
        "weight to obtain the absolute mL/h this column carries. The absolute cohort",
        "ultrafiltration rate in supplementary Table 3 is 1152 (773) mL/h. Meaningful only",
        "when RRT_CRRT_ACTIVE = 1; converted to L/h inside model(). Effluent flow was",
        "separately tested as a covariate ON the clearance arm and not retained -- Methods:",
        "'The surface area and the running time of the hemofilter were the covariates",
        "investigated for CLRRT', and only surface area survived. Founding example for this",
        "canonical column (operator ruling, sidecar request-001, answered 2026-08-28)."
      ),
      source_name        = "phi_Effluent"
    ),
    URINE_FLOW = list(
      description        = "Urine flow rate over the urine recovery interval containing the current observation",
      units              = "mL/h",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Does two jobs, both keyed on the register's documented zero sentinel ('0 = no urine",
        "collected during the interval'). (1) It is the flow denominator of the urine",
        "concentration observable. The supplementary material states that 'CUr were obtained",
        "as AUr/UVol, where UVol is the measured volume of urine excreted in each urine",
        "recovery interval', with AUr the amount accumulated over that same interval by",
        "eq. 3, dAUr/dt = (CLRenal/V1) * AP,Pre. Dividing an interval-cumulative amount by",
        "an interval-cumulative volume is algebraically the ratio of the two rates, so the",
        "observable is encoded here in the equivalent rate form CUr = CLRenal * CPre /",
        "phi_Urine -- exactly parallel to eq. 5 for the effluent, and numerically identical",
        "to the paper's own expression whenever clearance and urine flow are constant across",
        "the recovery interval, which is what a single measured UVol per interval already",
        "assumes. The reuse of URINE_FLOW in this flow-denominator role rather than a new",
        "canonical is an operator ruling (sidecar request-001, question 2, answered",
        "2026-08-28). (2) The same zero sentinel supplies the paper's no-diuresis gate on",
        "the renal arm: 'CLRenal and CLRRT were absent in patients without diuresis or",
        "without hemofilter' (supplementary material). Only 3 of the 13 CKRT patients showed",
        "residual diuresis (Table 1 footnote b). Table 1 reports 24-hour urine output as",
        "mean (SD) 981 (616) mL without CKRT and 184 (216) mL with CKRT; divide by 24 to",
        "obtain the mL/h this column carries when a finer per-interval measurement is not",
        "available."
      ),
      source_name        = "UVol"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "months",
      type        = "continuous",
      notes       = "Screened and not retained (Butragueno-Laiseca 2022 Results, Population pharmacokinetic modelling: 'No statistical correlations between age, albumin plasma levels, hematocrit and other laboratory values, and the PK parameters were found (p > 0.05)'). Median (range) 7 months (3 months to 15 years) overall; mean (SD) 29 (43) months without CKRT and 62 (62.5) months with CKRT (Table 1). The Discussion attributes the null postmenstrual-age result to the cohort composition -- 'half of the population was 2 years or older, and the minimum age was three months. The absence of premature or neonate patients might have precluded finding a relationship between PMA and CL' -- and offers the retained HT covariate as the organ-maturation marker in its place."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/dL",
      type        = "continuous",
      notes       = "Screened and not retained, despite piperacillin's plasma-protein binding. Mean (SD) 3.6 (0.5) g/dL without CKRT and 3.3 (0.6) g/dL with CKRT (Table 1)."
    ),
    HCT = list(
      description = "Haematocrit",
      units       = "%",
      type        = "continuous",
      notes       = "Screened and not retained on any structural parameter. Unlike the sibling ButraguenoLaiseca_2025_teicoplanin model, haematocrit does not enter the post-filter plasma-flow correction here either, because this paper estimated the blood-to-plasma ratio directly and found it 'not significantly different from 1' (Results), leaving phi_Pl,corr equal to the raw blood flow. Mean (SD) 31.5 (5.4)% without CKRT and 30.8 (5.5)% with CKRT (Table 1)."
    ),
    FILT_RUNTIME = list(
      description = "Haemofilter running time at the current observation",
      units       = "h",
      type        = "continuous",
      notes       = "Screened on the haemofilter clearance arm and not retained; only surface area survived (Methods, Pharmacokinetic analysis: 'The surface area and the running time of the hemofilter were the covariates investigated for CLRRT'; Results retain only surface area). Supplementary Table 3 reports mean (SD) running times of 56.3 (37.6) h for the 0.2 m2 filter, 60.5 (52.7) h for 0.6 m2 and 10.5 (13.4) h for 1.2 m2. Not a registered canonical column, and none is proposed, because the model does not use it; the closest registered relative is the FILT_AGE_HI above-48-hour indicator."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 32L,
    n_studies      = 1L,
    n_observations = 429L,
    age_range      = "3 months to 15 years",
    age_median     = "7 months",
    weight_range   = "4 to 63 kg",
    weight_median  = "8.1 kg",
    race_ethnicity = "Not reported.",
    disease_state  = paste(
      "Critically ill children admitted to a paediatric intensive care unit and treated",
      "with piperacillin-tazobactam, 13 of them undergoing continuous kidney replacement",
      "therapy in continuous venovenous haemodiafiltration modality (10 of the 13 with",
      "citrate anticoagulation). CKRT indications were acute kidney injury with oliguria",
      "in 10 patients and acute kidney injury with hypervolaemia in 2. Most patients were",
      "cardiac: 84% of each group carried a diagnosed cardiopathy, and 42% (non-CKRT) and",
      "77% (CKRT) were in the post-operative period. ALL patients were mechanically",
      "ventilated and received vasopressors. PRISM III score median 5 (non-CKRT) and 8.5",
      "(CKRT); ECMO in 1/19 and 3/13; mortality 26% and 46%. Supplementary Tables 1 and 2",
      "list the per-patient main diagnosis and isolated organism; treatment was empiric in",
      "the majority, and the organisms isolated were Pseudomonas aeruginosa, Klebsiella",
      "pneumoniae, Klebsiella oxytoca, Enterobacter cloacae and Stenotrophomonas",
      "maltophilia, with MICs of 8 to 64 mg/L."
    ),
    dose_range     = paste(
      "All patients received 100 mg/kg of piperacillin/tazobactam intravenously every 8 h.",
      "In patients with CKRT the dosing interval was increased to 12 h at the fourth dose,",
      "following the recommended renal adjustment. Sampling began once patients had",
      "received at least three doses, with blood and urine collected before (T0) and 2, 4,",
      "6 and 8 h after the start of the infusion; in CKRT patients, pre-filter, post-filter",
      "and effluent-port samples were drawn simultaneously from the Prismaflex device.",
      "Approximately 5 to 10 plasma and 5 urine samples per patient. Sampling times",
      "differed across patients (supplementary Figures 1 and 2). The infusion duration is",
      "not stated in the Methods; the 30-minute infusion used throughout the",
      "model-based dosing evaluation (supplementary Table 4, 'Current' schedules) is the",
      "value adopted in the validation vignette."
    ),
    regions        = "Hospital General Universitario Gregorio Maranon, Madrid, Spain (single centre).",
    renal_function = paste(
      "Non-CKRT group: eGFR (Schwartz) median 119.3 mL/min/1.73 m2, mean (SD) 123 (55);",
      "mean serum creatinine 0.37 (SD 0.19) mg/dL; 24-hour urine output 981 (616) mL. CKRT",
      "group: eGFR not measured; 24-hour urine output 184 (216) mL, calculated from the",
      "three patients who retained residual diuresis. CKRT prescription (supplementary",
      "Table 3): blood flow 70 (31) mL/min, citrate anticoagulation in 10 of 13 at a dose",
      "of 2.5 (0.5) mmol/L with a citrate flow of 655 (345) mL/h, substitution flow 125",
      "(200) mL/h, ultrafiltration rate 1152 (773) mL/h, extraction rate 59 (37) mL/h, and",
      "heparin 20 (10) IU/kg/h in the remaining 3. Filter surface areas 0.2 m2 (n = 6),",
      "0.6 m2 (n = 4) and 1.2 m2 (n = 3), with effluent flows of 78 (24), 54 (5) and 60",
      "(15) mL/kg/h and mean running times of 56.3, 60.5 and 10.5 h respectively."
    ),
    notes          = paste(
      "Baseline demographics from Butragueno-Laiseca 2022 Table 1, reported separately for",
      "the 19 patients without and the 13 patients with CKRT; age, weight and height did",
      "not differ significantly between the groups. Sex is not reported. 429 piperacillin",
      "concentrations across five matrices (plasma, pre-filter, post-filter, effluent and",
      "urine) were modelled simultaneously. Concentrations were assayed by validated HPLC",
      "with UV detection at 254 nm, linear over 0.5 to 1000 mg/L with a limit of",
      "quantification of 0.5 mg/L (supplementary material, Analytical method); urine",
      "samples were diluted 1/20 or 1/50 before precipitation. Tazobactam was not measured.",
      "NONMEM 7.4 with FOCE-I on log-transformed concentrations; covariates selected by",
      "stepwise covariate modelling (forward 0.05, backward 0.01) with continuous",
      "covariates normalised by the population median; parameter uncertainty from 500",
      "bootstrap datasets; model evaluation by prediction-corrected VPC over 1000 simulated",
      "datasets."
    )
  )

  ini({
    # -----------------------------------------------------------------------
    # STRUCTURAL PARAMETERS -- Butragueno-Laiseca 2022 Table 2 (selected
    # model), point estimates with the bootstrap 95% confidence interval in
    # parentheses. Total elimination from the central compartment is the
    # ADDITIVE sum of three arms (supplementary material eq. 1):
    #
    #   dA_P,Pre/dt = (CL_D/V2)*A_2 - (CL_D/V1)*A_P,Pre
    #                 - ((CL_Renal + CL_RRT + CL_M)/V1)*A_P,Pre
    #
    # This is the substantive structural difference from the sibling
    # ButraguenoLaiseca_2025_teicoplanin model, whose renal and haemofilter
    # arms are mutually EXCLUSIVE (teicoplanin CL_R was estimated to be zero
    # on CKRT). Here they coexist, and a non-renal arm is additionally
    # identified: "Data supported the characterization of nonrenal clearance
    # (CL M) as an additional elimination pathway (p < 0.01)" (Results).
    #
    # lcl_hemodialysis is the typical CKRT clearance for the LARGE (1.2 m^2)
    # haemofilter, which is this paper's reference level of theta_FILT.
    # -----------------------------------------------------------------------
    lcl_renal        <- log(1.3);  label("Typical renal clearance CLR at median eGFR and height (L/h)")            # Table 2: theta_CLR = 1.3 L/h (bootstrap 95% CI 0.9-1.7)
    lcl_nonren       <- log(0.5);  label("Typical non-renal clearance CLM at 8.1 kg (L/h)")                        # Table 2: theta_CLM = 0.5 L/h (0.3-0.8)
    lcl_hemodialysis <- log(1.34); label("Typical haemofilter clearance CLCKRT with a large 1.2 m2 filter (L/h)")   # Table 2: theta_CLCKRT = 1.34 L/h (1.2-1.5)
    lvc              <- log(3.23); label("Typical central volume of distribution V1 at 8.1 kg (L)")                 # Table 2: theta_V = 3.23 L (CI printed as "2.6-0.8"; see vignette Errata)
    lvp              <- log(1.45); label("Typical peripheral volume of distribution V2 at 8.1 kg (L)")              # Table 2: theta_VT = 1.45 L (0.9-2.7)
    lq               <- log(0.51); label("Typical inter-compartmental clearance CLD (L/h)")                         # Table 2: theta_CLD = 0.51 L/h (0.4-1.2)

    # -----------------------------------------------------------------------
    # HAEMOFILTER-SIZE EFFECT ON THE CKRT CLEARANCE ARM -- Table 2:
    # CLCKRT = theta_CLCKRT x theta_FILT, with theta_FILT taking one value per
    # filter surface area and the LARGE filter as the reference:
    #
    #   High   (1.2 m^2) = 1     Fixed  (reference)
    #   Med    (0.6 m^2) = 0.74  (0.6-0.9)
    #   Low    (0.2 m^2) = 0.28  (0.2-0.4)
    #
    # Note the reference level is the opposite of the sibling teicoplanin
    # model's. The registered FILT_SA_MED / FILT_SA_LARGE columns keep their
    # canonical semantics (both 0 => small filter) and the small-filter
    # indicator is derived inside model() as (1 - FILT_SA_MED - FILT_SA_LARGE),
    # so all three printed Table 2 numbers appear verbatim and no
    # reference-category flip is imposed on the covariate register.
    #
    # The reference multiplier itself is a structural 1 that the authors fixed
    # ("theta_FILT_High = 1 Fixed"), so it is not carried as an ini() entry.
    # -----------------------------------------------------------------------
    e_filt_sa_med_cl_hemodialysis   <- 0.74; label("Multiplicative factor on CLCKRT for a medium 0.6 m2 filter (vs large)")  # Table 2: theta_FILT_Med = 0.74
    e_filt_sa_small_cl_hemodialysis <- 0.28; label("Multiplicative factor on CLCKRT for a small 0.2 m2 filter (vs large)")   # Table 2: theta_FILT_Low = 0.28

    # -----------------------------------------------------------------------
    # INTERPATIENT VARIABILITY -- Table 2, "IPV" column.
    #
    # Methods: "Interpatient variability was modelled using the exponential
    # model". Table 2 footnote: "IPV is expressed as CV (coefficient of
    # variation, %) calculated as sqrt(exp(omega^2) - 1) x 100, where omega^2
    # corresponds to the variance of the random effects." Inverting that
    # relation gives omega^2 = log(1 + CV^2), the standard lognormal
    # conversion.
    #
    #   IPV CLR    =  70 % CV (bootstrap 95% CI 40-96),  eta-shrinkage 28%
    #   IPV CLCKRT =  20 % CV (5-26),                    eta-shrinkage 45%
    #   IPV CLM    = 123 % CV (59-239),                  eta-shrinkage 21%
    #   IPV V1     =  29 % CV (11-44),                   eta-shrinkage 27%
    #   IPV V2     =  91 % CV (28-188),                  eta-shrinkage 33%
    #
    # No IPV on CLD: Table 2 records "Ne" (not estimated) for both its IPV and
    # its shrinkage. The etas are independent -- the supplementary material
    # records that "the significance of the diagonal and off-diagonal elements
    # of the variance-covariance matrix" was explored, and Table 2 reports no
    # off-diagonal element.
    # -----------------------------------------------------------------------
    etalcl_renal        ~ log(1 + 0.70^2)  # Table 2: IPV CLR    =  70 % CV
    etalcl_hemodialysis ~ log(1 + 0.20^2)  # Table 2: IPV CLCKRT =  20 % CV
    etalcl_nonren       ~ log(1 + 1.23^2)  # Table 2: IPV CLM    = 123 % CV
    etalvc              ~ log(1 + 0.29^2)  # Table 2: IPV V1     =  29 % CV
    etalvp              ~ log(1 + 0.91^2)  # Table 2: IPV V2     =  91 % CV

    # -----------------------------------------------------------------------
    # RESIDUAL ERROR -- Table 2, "Residual error (Log (mg/L))" row.
    #
    # Methods: "Piperacillin concentrations were logarithmically transformed
    # for the analysis. ... the residual error was characterized with an
    # additive model on the logarithmic scale." Supplementary material:
    # "Different magnitudes of residual error were estimated for each of the
    # different types of measured entities."
    #
    # An additive residual on log-transformed concentrations is a lognormal
    # residual on the linear concentration scale, so all four magnitudes map to
    # the canonical expSd / lnorm() form rather than to propSd / addSd.
    #
    #   Plasma    0.28 (0.2-0.4),   eps-shrinkage 11%
    #   Postfilter 0.17 (0.1-0.2),  eps-shrinkage 12%
    #   Effluent  0.35 (0.15-0.5),  eps-shrinkage  6%
    #   Urine     1.2  (0.9-1.5),   eps-shrinkage 34%
    #
    # Table 2 footnote e attaches to the "Plasma" label and reads "Include
    # pre-filter concentrations measured in patients with continuous kidney
    # replacement therapy (CKRT)", so plasma and pre-filter share one residual
    # magnitude and are both represented by the Cc observable -- consistent
    # with the supplementary material, where the central compartment is named
    # A_P,Pre and "Predicted CP and CPre were obtained as A1/V1".
    # -----------------------------------------------------------------------
    expSd             <- 0.28; label("Log-scale additive residual SD for plasma and pre-filter concentrations (unitless)")  # Table 2: plasma residual error 0.28
    expSd_Cpostfilter <- 0.17; label("Log-scale additive residual SD for post-filter concentrations (unitless)")            # Table 2: postfilter residual error 0.17
    expSd_Ceffluent   <- 0.35; label("Log-scale additive residual SD for effluent concentrations (unitless)")               # Table 2: effluent residual error 0.35
    expSd_Curine      <- 1.2;  label("Log-scale additive residual SD for urine concentrations (unitless)")                  # Table 2: urine residual error 1.2
  })

  model({
    # -------------------------------------------------------------------
    # 1. Individual PK parameters.
    #
    # RENAL ARM. Table 2 gives CLR = theta_CLR * (eGFR/119.3) * (HGT/69),
    # both ratios entering LINEARLY -- there are no exponents anywhere in
    # Table 2. Two independent gates apply:
    #
    #   * The eGFR and height covariates apply only OFF CKRT. Table 2
    #     footnote a: "In CKRT patients showing residual diuresis, CLR is
    #     given solely by the typical estimate theta_CLR." eGFR was not
    #     measured in CKRT patients at all (Table 1 prints a dash).
    #     renal_cov collapses to exactly 1 when RRT_CRRT_ACTIVE = 1.
    #
    #   * The whole arm is gated off in a patient with no diuresis
    #     (supplementary material: "CLRenal and CLRRT were absent in
    #     patients without diuresis or without hemofilter"), carried by the
    #     URINE_FLOW zero sentinel that the covariate register already
    #     documents for that column.
    # -------------------------------------------------------------------
    renal_cov <- (1 - RRT_CRRT_ACTIVE) * (CRCL / 119.3) * (HT / 69) + RRT_CRRT_ACTIVE

    cl_renal <- (URINE_FLOW > 0) * exp(lcl_renal + etalcl_renal) * renal_cov

    # HAEMOFILTER ARM. Table 2: CLCKRT = theta_CLCKRT * theta_FILT, gated on
    # by RRT_CRRT_ACTIVE. The large filter is the reference (multiplier 1),
    # so only the medium and small multipliers appear; the small-filter
    # indicator is derived from the two registered canonical columns.
    filt_sa_small <- 1 - FILT_SA_MED - FILT_SA_LARGE

    cl_hemodialysis <-
      RRT_CRRT_ACTIVE * exp(lcl_hemodialysis + etalcl_hemodialysis) *
      e_filt_sa_med_cl_hemodialysis^FILT_SA_MED *
      e_filt_sa_small_cl_hemodialysis^filt_sa_small

    # NON-RENAL ARM and the volumes. All three carry the plain linear weight
    # ratio (WGT/8.1) with NO allometric exponent -- Results: "Adding the
    # allometric scaling factor to the CLM vs Body weight relationship did
    # not improve the fit (p > 0.05), and fixing it to 0.75 increased the
    # -2LL value." CLD carries no covariate and no IPV.
    cl_nonren <- exp(lcl_nonren + etalcl_nonren) * (WT / 8.1)

    cl <- cl_renal + cl_nonren + cl_hemodialysis
    vc <- exp(lvc + etalvc) * (WT / 8.1)
    vp <- exp(lvp + etalvp) * (WT / 8.1)
    q  <- exp(lq)

    # 2. Micro-constants.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # -------------------------------------------------------------------
    # 3. ODE system -- supplementary material eq. 1-3, verbatim. The central
    #    compartment is the authors' A_P,Pre: it is both the plasma and the
    #    pre-filter compartment. `urine` is their A_Ur, the cumulative
    #    amount renally excreted, driven by the renal arm alone.
    #
    #    Piperacillin was given as an intravenous infusion; the infusion
    #    duration is supplied on the dose records rather than modelled,
    #    since the paper does not estimate it.
    # -------------------------------------------------------------------
    d/dt(central)     <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <-   k12 * central - k21 * peripheral1
    d/dt(urine)       <-   cl_renal / vc * central

    # -------------------------------------------------------------------
    # 4. Observations -- four matrices, four residual magnitudes.
    #
    # Cc is the plasma concentration and, in CKRT patients, the pre-filter
    # concentration; the two are the same model prediction and share one
    # residual magnitude (Table 2 footnote e).
    #
    # POST-FILTER (supplementary eq. 4): CPost = CPre * (1 - CLRRT /
    # phi_Pl,corr), with phi_Pl,corr = phi_Blood x BPR. The Results pin
    # BPR: "The parameter blood to plasma ratio was not significantly
    # different from 1 (p < 0.01)", so phi_Pl,corr is just the blood flow
    # and no haematocrit correction is applied (contrast the sibling
    # teicoplanin model, where BPR was not estimated and the plasma
    # fraction 1 - HCT/100 stands in for it).
    #
    # BFR is stored in the canonical mL/min and converted to L/h (x 0.06).
    # That per-minute reading, not the main text's "5 (SD 2) ml/kg/h", is
    # what reproduces the paper's own published extraction ratio of 0.18
    # (range 0.11-0.25): at the supplementary Table 3 blood flow of
    # 70 mL/min = 4.2 L/h, the three filter sizes give CLRRT/phi_Blood of
    # 0.38/4.2 = 0.09, 1.01/4.2 = 0.24 and 1.34/4.2 = 0.32, whose
    # filter-count-weighted mean over the 6/4/3 patients is 0.19.
    #
    # EFFLUENT (supplementary eq. 5): CEffl = CLRRT * CPre / phi_Effl.
    # Together with the post-filter equation this reproduces the paper's
    # other published derived quantity, the sieving coefficient
    # SC = CEffl / ((CPre + CPost)/2), whose weighted mean over the same
    # patients is 0.81 against a published 0.77 (range 0.5-1).
    #
    # URINE: the supplementary material defines CUr = AUr/UVol over each
    # urine recovery interval. Dividing an interval-cumulative amount by an
    # interval-cumulative volume is the ratio of the two rates, so the
    # observable is written here in the equivalent rate form
    # CUr = CLRenal * CPre / phi_Urine, exactly parallel to eq. 5. The
    # `urine` state above still integrates eq. 3 verbatim and is what a
    # mass-balance check should use.
    #
    # Every denominator carries a small numerical floor so that a subject
    # off the circuit (BFR = 0, effluent flow = 0) or without diuresis
    # (URINE_FLOW = 0) cannot divide by zero, and so that no sampled
    # combination of a large filter with a low blood flow can drive the
    # post-filter concentration negative. None of the floors binds anywhere
    # in the covariate ranges the paper reports.
    # -------------------------------------------------------------------
    plasma_flow_filter <- max(BFR * 0.06, 0.001)
    effluent_flow      <- max(RRT_CRRT_EFFLUENT_FLOW / 1000, 0.001)
    urine_flow         <- max(URINE_FLOW / 1000, 0.001)

    filter_extraction <- cl_hemodialysis / plasma_flow_filter

    Cc          <- central / vc
    Cpostfilter <- Cc * max(1 - filter_extraction, 0.001)
    Ceffluent   <- cl_hemodialysis * Cc / effluent_flow
    Curine      <- cl_renal * Cc / urine_flow

    Cc          ~ lnorm(expSd)
    Cpostfilter ~ lnorm(expSd_Cpostfilter)
    Ceffluent   ~ lnorm(expSd_Ceffluent)
    Curine      ~ lnorm(expSd_Curine)
  })
}
