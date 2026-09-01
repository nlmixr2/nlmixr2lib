ButraguenoLaiseca_2024_meropenem <- function() {
  description <- paste(
    "Two-compartment population PK model for intravenous meropenem in 16",
    "critically ill children in a paediatric ICU, 7 of whom were receiving",
    "continuous kidney replacement therapy (CKRT) (Butragueno-Laiseca 2024).",
    "Five matrices were fitted simultaneously following the Broeker et al.",
    "approach -- plasma, pre-filter, post-filter, effluent and urine -- which",
    "untangles total elimination into three additive arms: a renal arm",
    "(CLR, 5.5 L/h at full maturity and a median eGFR of 95 mL/min/1.73 m2,",
    "encoded as lcl_renal), a non-renal / metabolic arm (CLM, 0.91 L/h,",
    "encoded as lcl_nonren) and a CKRT arm (CLCKRT, 1.26 L/h for the medium",
    "hemofilter, encoded as lcl_hemodialysis). The renal arm scales linearly",
    "with eGFR and carries a postnatal-age maturation function reaching half",
    "of the mature value at 15.6 months (Hill coefficient not significantly",
    "different from 1); in CKRT patients it is replaced by a flat 0.96 L/h",
    "for those retaining residual diuresis and switched off entirely in",
    "anuric patients. Hemofilter surface area is a categorical covariate on",
    "the CKRT arm (medium 0.6 m2 is the reference; low 0.2 m2 multiplier",
    "0.27, high 1.2 m2 multiplier 1.5). Central volume scales linearly with",
    "body weight; the peripheral volume and the inter-compartmental clearance",
    "carry no covariates and no IIV. Post-filter and effluent concentrations",
    "are derived algebraically from the pre-filter (central) concentration",
    "using the measured circuit blood and effluent flows, and the urinary",
    "concentration from a cumulative urinary-amount state divided by the",
    "collected volume; each of the four observables carries its own",
    "log-scale residual error."
  )
  reference <- paste(
    "Butragueno-Laiseca L, Troconiz IF, Grau S, Campillo N, Padilla B,",
    "Fernandez SN, Slocker M, Herrera L, Santiago MJ.",
    "How to use meropenem in pediatric patients undergoing CKRT? Integrated",
    "meropenem pharmacokinetic model for critically ill children.",
    "Antimicrob Agents Chemother. 2024;68(6):e0172923.",
    "doi:10.1128/aac.01729-23"
  )
  vignette <- "ButraguenoLaiseca_2024_meropenem"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central     = list(analyte = "meropenem", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "meropenem", units = "mg", specimen = "plasma", verified = TRUE),
    urine       = list(analyte = "meropenem", units = "mg", specimen = "urine",  verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "The only covariate retained on the central volume, and the only structural role body",
        "weight plays in this model: Butragueno-Laiseca 2024 Table 4 prints the parameter model",
        "for V1 as theta_V1 * WGT/10, a LINEAR (exponent 1) normalisation to the cohort median",
        "weight of 10 kg, which the same table's footnote states explicitly ('WGT, body weight",
        "(median = 10 kg)'). Results, Selection of covariates: 'Body weight shows significant",
        "covariate effects on V1 (P < 0.01) but did not affect any of the other model",
        "parameters (P > 0.05)', and the Discussion adds that weight was NOT selected on",
        "clearance in favour of postnatal age, 'reflecting organ maturation beyond the value of",
        "eGFR'. Median (IQR) weight was 7.5 (5.5-17.5) kg in the group without CKRT and",
        "20 (7.4-40) kg in the CKRT group (Table 1). Body weight is nevertheless strongly",
        "collinear with hemofilter size in this cohort because filters were assigned by weight",
        "band (low 3-10 kg, medium 10-30 kg, high 30-60 kg; Methods, Probability of target",
        "attainment), and with the effluent flow -- the Discussion attributes the difference",
        "from the Thy et al. model to exactly that correlation. Treated as a time-fixed",
        "baseline value in the source analysis."
      ),
      source_name        = "WGT"
    ),
    PNA = list(
      description        = "Postnatal (chronological) age",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Drives the renal maturation function on CLR in patients WITHOUT CKRT:",
        "Butragueno-Laiseca 2024 Table 4 prints the parameter model as",
        "theta_CLR * (eGFR/95) * AGE/(AGE + theta_50), an Emax (hyperbolic) form whose Hill",
        "coefficient the Results report as 'not significantly different from 1 (P > 0.05)'.",
        "Mapped to the canonical PNA rather than to AGE because the source is explicit that the",
        "quantity is POSTNATAL age standing in for organ maturation -- Discussion: 'weight was",
        "not selected in the final model in favor of post-natal age, reflecting organ maturation",
        "beyond the value of eGFR', and '90% of adult kidney functionality is predicted at ages",
        "over 12 months'. The canonical PNA is already in months, which is the unit the paper's",
        "theta_50 = 15.6 months is expressed in, so no conversion is needed inside model().",
        "Neonates under 28 days of postnatal age were excluded from the study (Limitations), so",
        "the maturation curve is extrapolated rather than fitted below about 1 month. Median",
        "(IQR) age was 8 (3.5-82) months without CKRT and 48 (5-106) months with CKRT",
        "(Table 1); the pooled cohort median implied by supplemental Table S1 is 19.5 months.",
        "The maturation term does NOT apply to CKRT patients, whose renal arm the paper reports",
        "as independent of age (Results, Selection of covariates)."
      ),
      source_name        = "AGE"
    ),
    CRCL = list(
      description        = "Estimated glomerular filtration rate (cystatin-C CKD-EPI and Schwartz equations), BSA-normalised",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Scales the renal clearance arm linearly through the ratio (CRCL/95), where 95",
        "mL/min/1.73 m^2 is the cohort median printed in the Butragueno-Laiseca 2024 Table 4",
        "footnote ('eGFR, estimated glomerular filtration rate (median = 95 mL/min/1.73 m2)').",
        "The linear, through-origin form is confirmed arithmetically by the Results: 'an",
        "increment in eGFR of 10 mL/min/1.73 m2 is associated with a 10.5% increase in CLR'",
        "-- 10/95 = 10.5% -- and again by the Discussion, 'the predicted percentage decreases",
        "in renal ... clearances for every reduction in 20 mL/min/1.73 m2 predicted by our",
        "model are 21%', since 20/95 = 21%. There is therefore no estimated coefficient on this",
        "covariate; the effect is the normalisation itself, so the reference value 95 appears as",
        "a hardcoded constant in model() rather than as an ini() parameter. Methods, Covariate",
        "selection: 'The eGFR was estimated using Cystatin C (Chronic Kidney Disease",
        "Epidemiology Collaboration [CKD-EPI] Equation) and Schwartz formula and integrated",
        "into the model for those patients without CKRT' -- accordingly the eGFR term is gated",
        "off for CKRT patients, whose Table 1 eGFR row reads 'N/A, not applicable'. Median (IQR)",
        "in the group without CKRT was 89 (69-122) mL/min/1.73 m^2."
      ),
      source_name        = "eGFR"
    ),
    RRT_CRRT_ACTIVE = list(
      description        = "CKRT-active indicator (1 while continuous kidney replacement therapy is running, 0 otherwise)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no CKRT running)",
      notes              = paste(
        "Selects between the two forms of the renal arm and gates the extracorporeal arm on.",
        "With RRT_CRRT_ACTIVE = 0 the renal arm takes its eGFR- and maturation-scaled form and",
        "CLCKRT is zero -- Methods, Base population model: 'CLCRT was absent in patients without",
        "hemofilters.' With RRT_CRRT_ACTIVE = 1 the renal arm collapses to the flat",
        "theta_CLR_CKRT = 0.96 L/h that 'resulted independently from eGFR and AGE (P > 0.05)'",
        "(Results, Selection of covariates) and the filter-size-adjusted CLCKRT switches on.",
        "The ACTIVE (time-varying) rather than the STATUS (subject-level) member of the",
        "RRT_<modality>_<kind> family is used for consistency with the sibling",
        "ButraguenoLaiseca_2025_teicoplanin.R from the same group and PICU, and because the",
        "quantity the model needs is whether the circuit is running at the current record; in",
        "this particular cohort every CKRT patient was on therapy throughout sampling, so the",
        "column does not in fact vary within subject in the source data. The modality was",
        "continuous kidney replacement therapy delivered on Prismaflex (Baxter Int.) devices",
        "with polyacrylonitrile AN69 hollow-fiber hemofilters (Results, Demographics and",
        "clinical data), hence CRRT rather than the intermittent-hemodialysis member."
      ),
      source_name        = "patient type (with or without CKRT)"
    ),
    URINE_VOL_24H = list(
      description        = "24-hour residual diuresis (total urine volume per day)",
      units              = "mL/24h",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Used ONLY as a binary preserved-diuresis gate on the CKRT renal arm, following the",
        "encoding ratified by Huppe_2023_fosfomycin.R. Butragueno-Laiseca 2024 reports the",
        "theta_CLR_CKRT = 0.96 L/h estimate specifically 'For the three patients under CKRT with",
        "diuresis' (Results, Selection of covariates), restated in the Discussion as 'for those",
        "patients that show diuresis, the estimate of renal clearance was five to sixfold lower",
        "than the corresponding value obtained in case of no CKRT'. The remaining four CKRT",
        "patients were anuric -- 'Four patients under CKRT had anuria, and the other three had",
        "urine volumes ranging from 10 to 123 mL during the time interval of the",
        "pharmacokinetic curve' (Results, Brief description of the data) -- and an anuric",
        "patient can have no urinary excretion, so the renal arm is switched off for them",
        "entirely. The gate is therefore written as (URINE_VOL_24H > 0), which is the source's",
        "own dichotomy (anuria versus any diuresis); note this differs from the 100 mL/24h",
        "clinical anuria cutoff Huppe 2023 uses, because this paper defines its groups by",
        "presence or absence of urine rather than by a volume threshold. Median (IQR) urine",
        "output was 141 (52-457) mL in the CKRT group and 1,100 (719-1,150) mL in the group",
        "without CKRT (Table 1). See the vignette Errata: the paper does not state that this",
        "gate is a model term, and the reading that the 0.96 L/h arm applies only to",
        "diuretic CKRT patients is an interpretation of the two quoted sentences."
      ),
      source_name        = "residual diuresis / urine output"
    ),
    URINE_VOL_INTERVAL = list(
      description        = "Urine volume collected in the urine recovery interval containing the current urinary observation",
      units              = "mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Converts the modelled cumulative urinary AMOUNT into the urinary CONCENTRATION the",
        "assay reports. Butragueno-Laiseca 2024 Methods, Base population model: 'CUr was",
        "obtained as AUr/UVol, where UVol is the measured volume of urine excreted in each urine",
        "recovery interval.' Urine was collected before (T0) and 2, 4, 6 and 8 h after the start",
        "of the infusion (Methods, Dosing, sampling, and analytical method), so a recovery",
        "interval is roughly 2 h and the urine state must be reset to zero at each interval",
        "boundary -- the validation vignette does this with evid = 6 records on the urine",
        "compartment. Stored in mL and divided by 1000 inside model() to give litres, so that the",
        "mg amount in the urine state yields mg/L. Not a covariate effect on any structural",
        "parameter: it appears only in the observation equation. Zero for an anuric interval, in",
        "which case there is no urinary observation to make; the division is floored in model()",
        "so that a zero value cannot produce a non-finite prediction. The paper reports urine",
        "output in patients without CKRT ranging from 14 to 125 mL/h, and 10-123 mL over the",
        "whole PK curve in the three CKRT patients with diuresis (Results, Brief description of",
        "the data), which brackets the plausible per-interval volumes. New canonical ratified by",
        "operator sidecar request-001 / response-001 question q2 option A."
      ),
      source_name        = "UVol"
    ),
    FILT_SA_MED = list(
      description        = "Medium haemofilter indicator (0.6 m2 membrane surface area)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (the low 0.2 m2 filter, the high 1.2 m2 filter, or no CKRT)",
      notes              = paste(
        "One of the three levels of the hemofilter-surface-area covariate on the CKRT clearance",
        "arm. IMPORTANT -- unlike the sibling ButraguenoLaiseca_2025_teicoplanin.R, this paper",
        "references the MEDIUM filter, not the small one: Butragueno-Laiseca 2024 Table 4 gives",
        "CLCKRT = theta_CLCKRT * theta_FILT with theta_FILT_Med = 1 (reference),",
        "theta_FILT_Low = 0.27 and theta_FILT_High = 1.5, and the Table 4 footnote reads",
        "'FILT_Low, Med, High filter surface 0.2, 0.6 (reference), and 1.2 m2, respectively'.",
        "The DATA columns are unchanged by that choice -- FILT_SA_MED and FILT_SA_LARGE keep",
        "their registered meaning of 'this filter is the medium / large one', with 0/0 denoting",
        "the small (low) filter -- and model() derives the low-filter indicator as",
        "1 - FILT_SA_MED - FILT_SA_LARGE so the published thetas can be carried verbatim",
        "(operator sidecar request-001 / response-001 question q3 option A). The multiplier",
        "attached to FILT_SA_MED is therefore the reference 1 and does not appear as an ini()",
        "parameter; the covariate is still referenced in model() through the derived low-filter",
        "indicator. The absolute sizes match the sibling paper exactly: HF20 0.2 m^2 in children",
        "under 10 kg, M60 0.6 m^2 between 10 and 30 kg, and M100 over 30 kg. Note the internal",
        "inconsistency in the source over the largest filter's area, recorded in the vignette",
        "Errata: Table 3 and the Table 4 footnote both say 1.2 m^2 while the Results text",
        "describes the M100 as 0.9 m^2. Two of the seven CKRT patients used the high filter,",
        "two the medium and three the low (Table 3). Meaningful only when RRT_CRRT_ACTIVE = 1."
      ),
      source_name        = "FILT (Med)"
    ),
    FILT_SA_LARGE = list(
      description        = "Large (high) haemofilter indicator (1.2 m2 membrane surface area)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (the medium 0.6 m2 reference filter, the low 0.2 m2 filter, or no CKRT)",
      notes              = paste(
        "The high-surface-area level of the hemofilter covariate (Butragueno-Laiseca 2024",
        "Table 4: theta_FILT_High = 1.5, RSE 23%), corroborated by Results, Selection of",
        "covariates: 'The depurative efficiency of the extracorporeal elimination was reduced",
        "and augmented by 73% and 50%, respectively, for low and high surface area filters",
        "compared to middle size' -- 0.27 is a 73% reduction and 1.5 a 50% augmentation.",
        "Mutually exclusive with FILT_SA_MED; both indicators are 0 for the low 0.2 m^2 filter,",
        "which in THIS paper is a non-reference level carrying the 0.27 multiplier rather than",
        "the reference. See the FILT_SA_MED entry for the reference-flip rationale. Meaningful",
        "only when RRT_CRRT_ACTIVE = 1."
      ),
      source_name        = "FILT (High)"
    ),
    BFR = list(
      description        = "Blood flow rate through the CKRT extracorporeal circuit",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters ONLY the post-filter observation equation, not the clearance model.",
        "Butragueno-Laiseca 2024 Methods, Base population model, equation 4:",
        "CPost = CPre * (1 - CLRRT / phi_Pl,corr), where 'The value of phi_Pl,corr was",
        "calculated as phi_Blood x BPR, phi_Pl,corr being the blood flow rate and BPR the blood",
        "to plasma concentration ratio, an additional parameter to be estimated from the model.'",
        "Median (IQR) blood flow was 70 (40-100) mL/min across the seven CKRT patients",
        "(Table 2). Stored in the canonical mL/min and converted to L/h (x 0.06) inside model()",
        "to match the L/h in which CLCKRT is expressed. Meaningful only when",
        "RRT_CRRT_ACTIVE = 1."
      ),
      source_name        = "Blood flow"
    ),
    HCT = list(
      description        = "Haematocrit",
      units              = "%",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Converts the extracorporeal blood flow to the corrected plasma flow entering the",
        "hemofilter in the post-filter observation equation. Unlike the sibling teicoplanin",
        "paper, where the plasma-fraction reading had to be inferred, THIS paper resolves the",
        "blood-to-plasma ratio for us -- Results, Base population model: 'the distribution of",
        "meropenem into the red blood cells was negligible, and the value of BPR shows a value",
        "equal to 1 hematocrit', i.e. BPR = 1 - haematocrit, the plasma volume fraction of",
        "blood. (The minus sign is dropped by the publisher's symbol font; the sentence is only",
        "arithmetically meaningful with it, and negligible red-cell distribution is precisely",
        "the condition under which the blood-to-plasma ratio equals the plasma fraction.) BPR",
        "is therefore not a free parameter in this implementation but the data-derived quantity",
        "1 - HCT/100. Haematocrit is not tabulated for this cohort; supply a physiological",
        "paediatric-ICU value. See the vignette Errata."
      ),
      source_name        = "hematocrit"
    ),
    RRT_CRRT_EFFLUENT_FLOW = list(
      description        = "Total effluent flow of the CKRT circuit (dialysate plus substitution plus net ultrafiltration)",
      units              = "mL/h",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters ONLY the effluent observation equation. Butragueno-Laiseca 2024 Methods, Base",
        "population model, equation 5: CEffl = CLRRT * CPre / phi_Effl, 'where phi_Pl,corr and",
        "phi_Effl are the corrected plasma and total effluent flows ... The values of phi_Blood",
        "and phi_Effluent were measured during the course of the study.' Table 2 reports the",
        "total median ultrafiltration rate as 1,253 (IQR 335-1,870) mL/h, with the dialysis",
        "(400, IQR 150-1,000 mL/h) and substitution (50, IQR 50-140 mL/h) legs listed",
        "separately; Table 3 reports the mean effluent flow per filter size as 44 (SD 8),",
        "60 (SD 4) and 46 (SD 12) mL/kg/h for the low, medium and high filters. Stored in mL/h",
        "and divided by 1000 inside model() to give L/h. The Discussion stresses that this is a",
        "measured, controlled circuit setting rather than a fitted quantity: 'meropenem",
        "concentrations in the effluent depend on CKRT parameters and also on the total effluent",
        "flow rate (Qeff), which is a variable that was controlled and known during the course",
        "of the treatment.' It was NOT retained as a covariate on any structural parameter --",
        "'neither body weight nor Qeff was included as significant covariates of CL' -- the",
        "authors attributing that to its high correlation with body weight. Meaningful only when",
        "RRT_CRRT_ACTIVE = 1; the division is floored in model() so a zero value cannot produce",
        "a non-finite prediction. New canonical ratified by operator sidecar request-001 /",
        "response-001 question q1."
      ),
      source_name        = "phi_Effl / Qeff"
    )
  )

  covariatesDataExcluded <- list(
    HT = list(
      description = "Body height",
      units       = "cm",
      type        = "continuous",
      notes       = "Screened in the stepwise covariate model and not retained. Butragueno-Laiseca 2024 Methods, Covariate selection lists height among the covariates 'tested in all pharmacokinetic parameters associated with interindividual variability'; Results, Selection of covariates report that 'Total protein or albumin serum levels, other demographics, and laboratory indexes did not impact model parameters.' Not tabulated for this cohort."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened and not retained (Butragueno-Laiseca 2024 Methods, Covariate selection lists albumin concentration; Results, Selection of covariates: 'Total protein or albumin serum levels ... did not impact model parameters'). A null albumin effect is expected for meropenem, whose unbound fraction in plasma the paper takes as 0.98 (Methods, Probability of target attainment). Not tabulated for this cohort."
    ),
    TPRO = list(
      description = "Total serum protein",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened and not retained alongside albumin (Results, Selection of covariates). Not tabulated for this cohort."
    ),
    BSA = list(
      description = "Body surface area",
      units       = "m^2",
      type        = "continuous",
      notes       = "Reported at baseline -- median (IQR) 0.38 (0.3-0.73) m^2 without CKRT and 0.79 (0.36-1.1) m^2 with CKRT (Table 1) -- and used to BSA-normalise the eGFR, but not retained as a covariate on any PK parameter in its own right. Distinct from the retained hemofilter membrane surface area (FILT_SA_MED / FILT_SA_LARGE), which is a device property rather than a body size."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 16L,
    n_studies      = 1L,
    n_observations = 212L,
    age_range      = "3 to 173 months (supplemental Table S1); neonates under 28 days of postnatal age were excluded",
    age_median     = "8 months (IQR 3.5-82) without CKRT; 48 months (IQR 5-106) with CKRT",
    weight_range   = "not reported directly; IQR 5.5-17.5 kg without CKRT and 7.4-40 kg with CKRT",
    weight_median  = "7.5 kg without CKRT; 20 kg with CKRT; 10 kg pooled (Table 4 footnote)",
    sex_female_pct = "Not reported.",
    race_ethnicity = "Not reported.",
    disease_state  = paste(
      "Critically ill children admitted to a paediatric intensive care unit and treated with",
      "meropenem, 7 of them undergoing continuous kidney replacement therapy on Prismaflex",
      "(Baxter Int.) devices with polyacrylonitrile AN69 hollow-fiber hemofilters. Most",
      "patients (62.5%) were in the postoperative period of congenital cardiopathy; 8/9",
      "without CKRT and 4/7 with CKRT carried a cardiopathy diagnosis, and 2/9 and 3/7",
      "respectively were on ECMO. Illness severity by PRISM III score was higher in the CKRT",
      "group (median 15, IQR 9.5-20) than without CKRT (5.5, IQR 0-12.25; P = 0.026).",
      "Mortality was 22% without and 57% with CKRT. Indications and isolates are listed in",
      "supplemental Table S1 (bacteraemia, septic shock, pneumonia, urinary tract infection;",
      "several patients received empiric therapy with no organism isolated). No patient",
      "experienced neurotoxicity or nephrotoxicity."
    ),
    dose_range     = paste(
      "Hospital protocol: 40 mg/kg intravenously every 8 h. In patients with CKRT the dosing",
      "interval was increased to every 12 h at the fourth dose, following the recommended",
      "renal adjustment; no patient in the no-CKRT group required adjustment. Infusions were",
      "given over 30 min (Results, Probability of target attainment, which describes the",
      "current schedule as '40 mg/kg/12 h from day 2 in a short infusion of 30 min').",
      "Sampling began once patients had received at least three doses, with blood and urine",
      "drawn before (T0) and 2, 4, 6 and 8 h after the start of the infusion. In CKRT",
      "patients, pre-filter, post-filter and effluent samples were drawn simultaneously,",
      "giving 18 samples per CKRT patient; patients without CKRT contributed 6 blood and 6",
      "urine samples."
    ),
    regions        = "Hospital General Universitario Gregorio Maranon, Madrid, Spain (single centre).",
    renal_function = paste(
      "Group without CKRT: eGFR median 89 (IQR 69-122) mL/min/1.73 m^2, urine output median",
      "1,100 (719-1,150) mL, hourly urine output 14-125 mL/h. CKRT group: eGFR not applicable;",
      "residual diuresis median 141 (52-457) mL, with four of the seven patients anuric and",
      "the other three producing 10-123 mL over the PK sampling window. CKRT prescription",
      "(Table 2): blood flow 70 (40-100) mL/min, citrate anticoagulation in 4 patients and",
      "continuous heparin in 3, citrate dose 2.8 (2.1-3.6) mmol/L, citrate flow 900",
      "(680-1,125) mL/h, substitution flow 50 (50-140) mL/h, dialysis flow 400 (150-1,000)",
      "mL/h, total ultrafiltration rate 1,253 (335-1,870) mL/h, extraction rate 80 (45-80)",
      "mL/h. Hemofilters (Table 3): low 0.2 m^2 HF20 (n = 3, mean effluent flow 44 mL/kg/h),",
      "medium 0.6 m^2 M60 (n = 2, 60 mL/kg/h), high 1.2 m^2 M100 (n = 2, 46 mL/kg/h)."
    ),
    notes          = paste(
      "Baseline demographics from Butragueno-Laiseca 2024 Table 1, reported separately for the",
      "9 patients without and the 7 patients with CKRT; age, weight, body surface area and",
      "both stay durations did not differ significantly between the groups. 212 meropenem",
      "concentrations were modelled, 118 (60%) of them from CKRT patients, split 99 plasma /",
      "pre-filter, 38 post-filter, 26 effluent and 49 urine; none was below the 0.5 mg/L limit",
      "of quantification. Concentrations were assayed by a validated HPLC-UV method (Waters",
      "Alliance e2695 with W2489 UV detection at 254 nm, Symmetry C18 column) linear from 0.5",
      "to 1,000 mg/L with accuracy below 11% and precision below 10% CV. NONMEM 7.4 with",
      "FOCE-I on log-transformed concentrations, PsN 4.9 for the stepwise covariate model,",
      "prediction-corrected VPCs and sampling importance resampling; the one-compartment model",
      "fitted significantly worse (P < 0.001) and the three-compartment model gave no",
      "significant improvement with non-identifiable parameters."
    )
  )

  ini({
    # -----------------------------------------------------------------------
    # STRUCTURAL PARAMETERS -- Butragueno-Laiseca 2024 Table 4 (selected
    # model). The paper's central claim is that fitting five matrices at once
    # untangles total elimination from the central compartment into three
    # ADDITIVE arms (Methods, Base population model, equation 1):
    #
    #   dA_P,Pre/dt = CLD/V2 * A2 - CLD/V1 * A_P,Pre
    #                 - (CLRenal + CLRRT + CLM)/V1 * A_P,Pre
    #
    #   lcl_renal        <- log(theta_CLR)    = log(5.5)   (Table 4)
    #   lcl_nonren       <- log(theta_CLM)    = log(0.91)  (Table 4)
    #   lcl_hemodialysis <- log(theta_CLCKRT) = log(1.26)  (Table 4)
    #
    # lcl_renal is the FULLY MATURE renal clearance at the median eGFR of
    # 95 mL/min/1.73 m^2, i.e. the asymptote of the maturation function --
    # Results, Selection of covariates: "The value of CLR for a patient with
    # normal and fully mature renal function is 5.5 L/h."
    #
    # lcl_hemodialysis is the CKRT clearance with the MEDIUM (0.6 m^2)
    # hemofilter, which is this paper's reference filter level. Note that the
    # sibling ButraguenoLaiseca_2025_teicoplanin.R references the SMALL filter
    # instead; the reference flip is why this model needs the derived
    # low-filter indicator below (operator sidecar q3 option A).
    #
    # `_hemodialysis` rather than `_dialysis` is used for the extracorporeal
    # arm to match the sibling teicoplanin model from the same group and PICU;
    # both suffixes are registered members of the multi-component-CL family.
    # -----------------------------------------------------------------------
    lcl_renal        <- log(5.5);  label("Typical renal clearance CLR, fully mature, at eGFR 95 mL/min/1.73 m2 (L/h)")  # Table 4: theta_CLR = 5.5 L/h (RSE 13%)
    lcl_nonren       <- log(0.91); label("Typical non-renal (metabolic) clearance CLM (L/h)")                            # Table 4: theta_CLM = 0.91 L/h (RSE 25%; SIR 95% CI 0.6-1.3)
    lcl_hemodialysis <- log(1.26); label("Typical CKRT clearance CLCKRT with a medium 0.6 m2 hemofilter (L/h)")          # Table 4: theta_CLCKRT = 1.26 L/h (RSE 19%; SIR 95% CI 0.9-1.6)
    lvc              <- log(4.75); label("Typical central volume of distribution V1 at 10 kg (L)")                       # Table 4: theta_V = 4.75 L (RSE 17%; SIR 95% CI 3.6-5.9)
    lvp              <- log(10.7); label("Typical peripheral volume of distribution V2 (L)")                             # Table 4: theta_VT = 10.7 L (RSE 41%; SIR 95% CI 6.6-17.6)
    lq               <- log(0.28); label("Typical inter-compartmental (distribution) clearance CLD (L/h)")               # Table 4: theta_CLD = 0.28 L/h (RSE 25%; SIR 95% CI 0.22-0.39)

    # -----------------------------------------------------------------------
    # RENAL MATURATION -- Butragueno-Laiseca 2024 Table 4 parameter model for
    # CLR:
    #
    #   CLR = theta_CLR * (eGFR/95) * AGE/(AGE + theta_50)
    #
    # an Emax (hyperbolic) maturation in postnatal age whose Hill coefficient
    # "was not significantly different from 1 (P > 0.05)" (Results, Selection
    # of covariates), hence no lhill parameter. theta_50 is log-transformed
    # here to keep it positive, following the Su_2016_dexmedetomidine.R
    # `tm50_cl` precedent for an Emax-form age maturation on a clearance.
    #
    # The eGFR term needs no ini() parameter: it is a bare normalisation by
    # the cohort median of 95 mL/min/1.73 m^2 with no estimated coefficient,
    # as the Results confirm arithmetically ("an increment in eGFR of
    # 10 mL/min/1.73 m2 is associated with a 10.5% increase in CLR"; 10/95 =
    # 10.5%). The constant 95 therefore appears in model() with a source
    # comment rather than as a parameter.
    #
    # ERRATUM: Table 4 prints theta_50 = 15.6 months while the Results text
    # of the same paper says "50% of adult renal function is achieved at the
    # age of 15.9 months". The Table 4 estimate is used, being the parameter
    # table; see the vignette Errata.
    # -----------------------------------------------------------------------
    tm50_cl_renal <- log(15.6); label("Postnatal age at which CLR reaches 50% of its mature value (months)")  # Table 4: theta_50 = 15.6 months (RSE 16%; SIR 95% CI 11.3-20)

    # -----------------------------------------------------------------------
    # RENAL CLEARANCE IN CKRT PATIENTS -- Table 4 lists a SECOND typical value
    # for the renal arm, applying to patients on CKRT:
    #
    #   CLR_CKRT (L/h)   theta_CLR_CKRT = 0.96 (53) (0.5-1.8)
    #
    # It is flat: Results, Selection of covariates -- "For the three patients
    # under CKRT with diuresis, CLR was 0.96 L/h (P < 0.01) and resulted
    # independently from eGFR and AGE (P > 0.05)." Encoded as a covariate
    # effect of CKRT status on the renal arm, carrying the published value in
    # absolute L/h (the same absolute-valued covariate-effect form the
    # register documents for Delattre 2010 amikacin) so that Table 4 is
    # transcribed verbatim rather than as a derived ratio.
    #
    # It applies only to CKRT patients who retain residual diuresis; the four
    # anuric patients of the seven have no renal elimination pathway at all
    # (see covariateData[[URINE_VOL_24H]] and the vignette Errata).
    # -----------------------------------------------------------------------
    e_rrt_crrt_active_cl_renal <- 0.96; label("Renal clearance CLR of a CKRT patient with residual diuresis (L/h)")  # Table 4: theta_CLR_CKRT = 0.96 L/h (RSE 53%; SIR 95% CI 0.5-1.8)

    # -----------------------------------------------------------------------
    # BODY-WEIGHT EFFECT ON THE CENTRAL VOLUME -- Table 4 prints the parameter
    # model for V1 as theta_V1 * WGT/10, i.e. a LINEAR normalisation to the
    # 10 kg cohort median (Table 4 footnote). The exponent of 1 is structural,
    # not estimated -- Table 4 reports no RSE or confidence interval for it --
    # so it is encoded with fixed(), matching the sibling teicoplanin model.
    # -----------------------------------------------------------------------
    e_wt_vc <- fixed(1); label("Exponent on (WT/10) for V1 (unitless)")  # Table 4 parameter model: theta_V1 * WGT/10, i.e. exponent 1

    # -----------------------------------------------------------------------
    # HEMOFILTER-SIZE EFFECT ON THE CKRT CLEARANCE ARM -- Table 4:
    # CLCKRT = theta_CLCKRT * theta_FILT, with theta_FILT taking one value per
    # filter surface area:
    #
    #   Low    (0.2 m^2) = 0.27  (RSE 22%)
    #   Medium (0.6 m^2) = 1     (reference)
    #   High   (1.2 m^2) = 1.5   (RSE 23%)
    #
    # corroborated by Results, Selection of covariates: "The depurative
    # efficiency of the extracorporeal elimination was reduced and augmented
    # by 73% and 50%, respectively, for low and high surface area filters
    # compared to middle size."
    #
    # The registered FILT_SA_MED / FILT_SA_LARGE data columns keep their
    # documented meaning (0/0 = the small/low filter), so model() derives the
    # low-filter indicator as 1 - FILT_SA_MED - FILT_SA_LARGE and applies the
    # two non-reference multipliers by the same power encoding the sibling
    # model uses. The medium multiplier is the reference 1 and needs no
    # parameter.
    # -----------------------------------------------------------------------
    e_filt_sa_small_cl_hemodialysis <- 0.27; label("Multiplicative factor on CLCKRT for a low 0.2 m2 filter (vs the medium reference)")   # Table 4: theta_FILT_Low = 0.27
    e_filt_sa_large_cl_hemodialysis <- 1.5;  label("Multiplicative factor on CLCKRT for a high 1.2 m2 filter (vs the medium reference)")  # Table 4: theta_FILT_High = 1.5

    # -----------------------------------------------------------------------
    # INTER-INDIVIDUAL VARIABILITY -- Table 4, "IPV (RSE%)" column.
    #
    # Methods: "Interpatient variability (IPV) was modeled using the
    # exponential model". Table 4 footnote: "IPV, inter-patient variability
    # expressed as CV (%) calculated as sqrt(e^(omega^2) - 1) x 100, where
    # omega^2 is the variance." Inverting gives omega^2 = log(1 + CV^2), the
    # standard lognormal conversion.
    #
    #   IPV CLR = 30 % CV (RSE 60%; SIR 95% CI 12-51), eta-shrinkage 43%
    #   IPV CLM = 86 % CV (RSE 59%; SIR 95% CI 60-122), eta-shrinkage 18%
    #   IPV V1  = 48 % CV (RSE 30%; SIR 95% CI 33-62),  eta-shrinkage 14%
    #
    # These are the FINAL (covariate) model values; Results, Selection of
    # covariates state that the covariates "reduced the estimates of the
    # inter-individual variability obtained from the base population to 30%,
    # 86%, and 48% for CLR, CLM, and V1, respectively" (from base-model
    # values of 115%, 81% and 103%).
    #
    # Table 4 reports IPV as "NS" (not significant) for CLCKRT, CLD and V2, so
    # those three carry no eta. The CLCKRT eta specifically vanished once
    # filter size was in the model -- Discussion: "once the surface area was
    # incorporated in the model as a covariate, the interindividual
    # variability in CLCKRT vanished." Covariances were not significant
    # (Results, Base population model), so the etas are independent.
    #
    # The single CLR eta is applied to the renal arm in BOTH its states: the
    # eGFR/maturation form and the flat CKRT form are two typical values of
    # one parameter, and Table 4 reports one IPV for the pair.
    # -----------------------------------------------------------------------
    etalcl_renal  ~ log(1 + 0.30^2)  # Table 4: IPV CLR = 30 % CV
    etalcl_nonren ~ log(1 + 0.86^2)  # Table 4: IPV CLM = 86 % CV
    etalvc        ~ log(1 + 0.48^2)  # Table 4: IPV V1  = 48 % CV

    # -----------------------------------------------------------------------
    # RESIDUAL ERROR -- Table 4, "Residual error [Log (mg/L)]" rows.
    #
    # Methods, Pharmacokinetic analysis: "Meropenem concentration data were
    # logarithmically transformed for the analysis. ... the residual error was
    # characterized with an additive model on the logarithmic scale. Different
    # magnitudes of residual error were estimated for each of the different
    # types of measured concentrations."
    #
    # An additive residual on log-transformed concentrations is a lognormal
    # residual on the linear scale, so all four map to the canonical expSd /
    # lnorm() form rather than to propSd / addSd.
    #
    #   Plasma and prefilter : 0.38 (RSE 12%; SIR 95% CI 0.34-0.45), eps-shrinkage 9%
    #   Postfilter           : 0.38 (RSE 30%; SIR 95% CI 0.31-0.52), eps-shrinkage 6%
    #   Effluent             : 0.37 (RSE 36%; SIR 95% CI 0.27-0.53), eps-shrinkage 4%
    #   Urine                : 0.86 (RSE 13%; SIR 95% CI 0.7-1.1),   eps-shrinkage 3%
    #
    # Plasma and pre-filter samples share one magnitude and are both
    # represented by the Cc observable, exactly as in the sibling model.
    # -----------------------------------------------------------------------
    expSd             <- 0.38; label("Log-scale additive residual SD for plasma and pre-filter concentrations (unitless)")  # Table 4: plasma and prefilter residual error 0.38
    expSd_Cpostfilter <- 0.38; label("Log-scale additive residual SD for post-filter concentrations (unitless)")            # Table 4: postfilter residual error 0.38
    expSd_Ceffluent   <- 0.37; label("Log-scale additive residual SD for effluent concentrations (unitless)")               # Table 4: effluent residual error 0.37
    expSd_Curine      <- 0.86; label("Log-scale additive residual SD for urine concentrations (unitless)")                  # Table 4: urine residual error 0.86
  })

  model({
    # -------------------------------------------------------------------
    # 1. Individual PK parameters.
    #
    # The renal arm has two mutually exclusive typical values selected by
    # CKRT status (Table 4 rows CLR and CLR_CKRT):
    #
    #   off CKRT : theta_CLR * (eGFR/95) * AGE/(AGE + theta_50)
    #   on  CKRT : theta_CLR_CKRT, flat, and only for patients who retain
    #              residual diuresis -- an anuric patient excretes nothing.
    #
    # 95 mL/min/1.73 m^2 is the cohort median eGFR printed in the Table 4
    # footnote; the normalisation carries no estimated coefficient.
    # -------------------------------------------------------------------
    maturation_cl_renal <- PNA / (PNA + exp(tm50_cl_renal))
    diuresis            <- (URINE_VOL_24H > 0)

    cl_renal_typ <-
      (1 - RRT_CRRT_ACTIVE) * exp(lcl_renal) * (CRCL / 95) * maturation_cl_renal +
      RRT_CRRT_ACTIVE * diuresis * e_rrt_crrt_active_cl_renal
    cl_renal <- cl_renal_typ * exp(etalcl_renal)

    cl_nonren <- exp(lcl_nonren + etalcl_nonren)

    # The low (0.2 m^2) filter is this paper's non-reference level but the
    # registered 0/0 state of the FILT_SA_* data columns, so it is derived
    # rather than supplied. The medium (0.6 m^2) filter is the reference and
    # contributes a multiplier of 1 through both exponents being 0.
    filt_sa_small <- 1 - FILT_SA_MED - FILT_SA_LARGE

    cl_hemodialysis <-
      RRT_CRRT_ACTIVE * exp(lcl_hemodialysis) *
      e_filt_sa_small_cl_hemodialysis^filt_sa_small *
      e_filt_sa_large_cl_hemodialysis^FILT_SA_LARGE

    cl <- cl_renal + cl_nonren + cl_hemodialysis
    vc <- exp(lvc + etalvc) * (WT / 10)^e_wt_vc
    vp <- exp(lvp)
    q  <- exp(lq)

    # 2. Micro-constants.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # 3. Concentration in the central compartment, which is also the
    #    pre-filter concentration (Methods, Base population model:
    #    "Predicted CP and CPre were obtained as A1/V1").
    Cc <- central / vc

    # -------------------------------------------------------------------
    # 4. ODE system -- Butragueno-Laiseca 2024 equations 1-3. Meropenem was
    #    given as a 30-min intravenous infusion; the infusion duration is
    #    supplied on the dose records rather than modelled, since the paper
    #    does not estimate it.
    #
    #    The urine state accumulates the AMOUNT excreted renally,
    #    dAUr/dt = CLRenal/V1 * A_P,Pre = cl_renal * Cc (equation 3). It must
    #    be reset to zero at each urine recovery-interval boundary; the
    #    validation vignette does this with `evid = 5, amt = 0, cmt = "urine"`
    #    replacement records placed a numerical epsilon AFTER each boundary,
    #    since rxode2 applies dose-type records before same-time observations
    #    and a reset sitting exactly on the boundary would zero the state
    #    before the end-of-interval sample is taken.
    # -------------------------------------------------------------------
    d/dt(central)     <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <-   k12 * central - k21 * peripheral1
    d/dt(urine)       <-   cl_renal * Cc

    # -------------------------------------------------------------------
    # 5. Observations.
    #
    # Post-filter (equation 4):
    #
    #   CPost = CPre * (1 - CLRRT / phi_Pl,corr),  phi_Pl,corr = phi_Blood * BPR
    #
    # The paper resolves BPR for us -- Results, Base population model: "the
    # distribution of meropenem into the red blood cells was negligible, and
    # the value of BPR shows a value equal to 1 [-] hematocrit" -- so the
    # corrected plasma flow is the blood flow times the plasma volume
    # fraction, BFR * (1 - HCT/100). BFR is stored in the canonical mL/min
    # and converted to L/h (x 0.06) to match the L/h of CLCKRT.
    #
    # Effluent (equation 5):
    #
    #   CEffl = CLRRT * CPre / phi_Effl
    #
    # with the total effluent flow stored in mL/h and converted to L/h.
    #
    # Urine: CUr = AUr / UVol, with the collected volume stored in mL and
    # converted to L so that the mg amount yields mg/L.
    #
    # Each denominator carries a small numerical floor so that non-CKRT
    # subjects (BFR = 0, effluent flow = 0), anuric intervals (collected
    # volume = 0) and any extreme sampled combination of a high-surface
    # filter with a low blood flow cannot produce a division by zero or a
    # negative concentration. None of the floors binds anywhere in the
    # covariate ranges the paper reports.
    # -------------------------------------------------------------------
    plasma_flow_filter <- max(BFR * 0.06 * (1 - HCT / 100), 0.001)
    effluent_flow      <- max(RRT_CRRT_EFFLUENT_FLOW / 1000, 0.001)
    urine_volume       <- max(URINE_VOL_INTERVAL / 1000, 0.001)

    Cpostfilter <- Cc * max(1 - cl_hemodialysis / plasma_flow_filter, 0.001)
    Ceffluent   <- cl_hemodialysis * Cc / effluent_flow
    Curine      <- urine / urine_volume

    Cc          ~ lnorm(expSd)
    Cpostfilter ~ lnorm(expSd_Cpostfilter)
    Ceffluent   ~ lnorm(expSd_Ceffluent)
    Curine      ~ lnorm(expSd_Curine)
  })
}
