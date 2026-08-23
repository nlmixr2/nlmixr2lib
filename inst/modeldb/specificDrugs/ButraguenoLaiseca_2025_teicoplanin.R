ButraguenoLaiseca_2025_teicoplanin <- function() {
  description <- paste(
    "Two-compartment population PK model for intravenous teicoplanin in 26",
    "critically ill children in a paediatric ICU, 12 of whom were receiving",
    "continuous kidney replacement therapy (CKRT) in continuous venovenous",
    "haemodiafiltration modality (Butragueno-Laiseca 2025). Plasma,",
    "pre-filter and post-filter concentrations were fitted simultaneously.",
    "Elimination from the central compartment is partitioned into a renal",
    "arm (CLR, 0.169 L/h typical at 8 kg, allometrically scaled with a fixed",
    "exponent of 0.75, encoded as lcl_renal) and a haemofilter arm (CLKRT,",
    "0.119 L/h typical for the small filter, encoded as lcl_hemodialysis).",
    "The two arms are mutually exclusive: renal clearance was estimated to",
    "be zero while CKRT was running, so total clearance is CLR off CKRT and",
    "CLKRT on CKRT. Haemofilter surface area is a categorical covariate on",
    "the CKRT arm (small 0.2 m2 reference; medium 0.6 m2 multiplier 3.58;",
    "large 1.2 m2 multiplier 5.04). Central volume scales allometrically",
    "with body weight with a fixed exponent of 1; the peripheral volume and",
    "the inter-compartmental clearance carry no covariates and no IIV. The",
    "post-filter concentration is derived algebraically from the pre-filter",
    "(central) concentration and the filter extraction ratio CLKRT divided",
    "by the plasma flow entering the haemofilter, and carries its own",
    "residual error. No other patient characteristic, including estimated",
    "glomerular filtration rate and serum albumin, was retained."
  )
  reference <- paste(
    "Butragueno-Laiseca L, Garcia-Orueta G, Riva N, Troconiz IF,",
    "Fernandez SN, Camacho Vicente V, Padilla B, Slocker M, Santiago MJ.",
    "Population pharmacokinetic analysis of teicoplanin in paediatric",
    "patients, including those receiving continuous kidney replacement",
    "therapy: a prospective cohort study.",
    "J Antimicrob Chemother. 2025;80(3):868-876.",
    "doi:10.1093/jac/dkaf012"
  )
  vignette <- "ButraguenoLaiseca_2025_teicoplanin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central     = list(analyte = "teicoplanin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "teicoplanin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric size descriptor on the renal clearance arm (fixed exponent 0.75) and on",
        "the central volume (fixed exponent 1), both normalised to a reference weight of 8 kg",
        "-- the cohort median (Butragueno-Laiseca 2025 Results, Patient population: median",
        "(range) weight 8 kg (4.3-44 kg); the (WGT/8) normalisation is printed in the",
        "Parameter model column of Table 2). The reference weight is confirmed arithmetically",
        "by the Discussion statement that the model predicts an increase of 0.016 L/h in CLR",
        "and 0.19 L in V1 per additional kg of body weight: 0.169 * 0.75 / 8 = 0.0158 L/h/kg",
        "and 1.56 / 8 = 0.195 L/kg. Body weight was NOT retained on the peripheral volume or",
        "on the inter-compartmental clearance (Results, Covariate selection). Treated as a",
        "time-fixed baseline value in the source analysis. Weight also acts as an implicit",
        "surrogate for haemofilter size, because filter size was assigned by weight band",
        "(small 3-10 kg, medium 10-30 kg, large 30-60 kg; Results, Patient population)."
      ),
      source_name        = "WGT"
    ),
    RRT_CRRT_ACTIVE = list(
      description        = "CKRT-active indicator (1 while continuous venovenous haemodiafiltration is running, 0 otherwise)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no CKRT running)",
      notes              = paste(
        "Switches the model between its two mutually exclusive elimination arms. Renal",
        "clearance is gated OFF while CKRT is running -- Table 2 footnote a: 'CLR in patients",
        "with CKRT was estimated to be zero', restated in the Supplementary material",
        "(Base population model): 'we have dismissed renal clearance (CLR) in patients with",
        "CKRT, as our analyzed data suggested that there was no renal elimination in these",
        "patients, thus, CL equals CLKRT'. The haemofilter arm is correspondingly gated ON.",
        "The ACTIVE (time-varying) rather than the STATUS (subject-level) member of the",
        "RRT_<modality>_<kind> family is used because the covariate genuinely varies within",
        "subject in the source data: Table 1 footnote a records that one patient started",
        "treatment without CKRT and required CKRT after the fifth dose, contributing to both",
        "groups. The modality is continuous venovenous haemodiafiltration (Results, Patient",
        "population), hence CRRT rather than the intermittent-haemodialysis member."
      ),
      source_name        = "CKRT"
    ),
    FILT_SA_MED = list(
      description        = "Medium haemofilter indicator (0.6 m2 membrane surface area)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (small 0.2 m2 filter, the reference category, or no CKRT)",
      notes              = paste(
        "One of the two non-reference levels of the three-level haemofilter-surface-area",
        "covariate on the CKRT clearance arm (Butragueno-Laiseca 2025 Table 2: CLKRT =",
        "theta_CLKRT * theta_FILT, with theta_FILT Small = 1 (reference), Medium = 3.58,",
        "Large = 5.04). Filter sizes were assigned by weight band -- small 3-10 kg, medium",
        "10-30 kg, large 30-60 kg (Results, Patient population) -- and Table 1 records 7",
        "small, 3 medium and 2 large filters among the 12 CKRT patients. Meaningful only when",
        "RRT_CRRT_ACTIVE = 1; the whole CKRT arm is gated off otherwise. Mutually exclusive",
        "with FILT_SA_LARGE. Effluent flow and blood flow were also tested on CLKRT but could",
        "not be retained alongside filter size because of their strong mutual correlation",
        "(Results, Covariate selection; Figure S3)."
      ),
      source_name        = "FILT (Med)"
    ),
    FILT_SA_LARGE = list(
      description        = "Large haemofilter indicator (1.2 m2 membrane surface area)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (small 0.2 m2 filter, the reference category, or no CKRT)",
      notes              = paste(
        "The second non-reference level of the haemofilter-surface-area covariate on the CKRT",
        "clearance arm (Butragueno-Laiseca 2025 Table 2: theta_FILT Large = 5.04, RSE 14%,",
        "95% CI 3.80-6.26). Meaningful only when RRT_CRRT_ACTIVE = 1. Mutually exclusive with",
        "FILT_SA_MED; both indicators are 0 for the small (0.2 m2) reference filter. Only 2 of",
        "the 12 CKRT patients used a large filter (Table 1)."
      ),
      source_name        = "FILT (Large)"
    ),
    BFR = list(
      description        = "Blood flow rate through the CKRT extracorporeal circuit",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters ONLY the post-filter observation equation, not the clearance model: the",
        "Supplementary material (Pharmacokinetic analysis) gives Cpost = Cpre * (1 - CLKRT /",
        "phi_PlCorr), where phi_PlCorr is 'the corrected plasma flow that goes into the",
        "hemofilter, which is calculated as phi_Blood x BPR (blood to plasma ratio)'. Blood",
        "flow was separately tested as a covariate ON CLKRT and rejected because it could not",
        "be retained alongside filter size (Results, Covariate selection). Median (range)",
        "50 (34-150) mL/min across the 12 CKRT patients (Table 1). Constant per subject",
        "during the study -- Supplementary material: 'For each individual, all the flow values",
        "(phi) remained constant during the study.' Meaningful only when RRT_CRRT_ACTIVE = 1;",
        "converted to L/h inside model()."
      ),
      source_name        = "Blood flow"
    ),
    HCT = list(
      description        = "Haematocrit",
      units              = "%",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Converts the extracorporeal blood flow to the plasma flow entering the haemofilter",
        "in the post-filter observation equation (plasma fraction of blood = 1 - HCT/100).",
        "Butragueno-Laiseca 2025 Results, Covariate selection state that no laboratory value",
        "significantly affected any PK parameter 'although CKRT parameters and haematocrit",
        "were included in the model'; haematocrit therefore enters the flow correction rather",
        "than any structural parameter. Median (range) 30.6% (20-38.7%) without haemofilter",
        "and 29.3% (26.3-37.5%) with haemofilter (Table 1). See the vignette Errata: the",
        "paper defines phi_PlCorr as phi_Blood x BPR but does not print the numerical",
        "definition of BPR; the plasma-fraction reading 1 - HCT/100 is the interpretation",
        "encoded here."
      ),
      source_name        = "Haematocrit"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "months",
      type        = "continuous",
      notes       = "Screened and not retained (Butragueno-Laiseca 2025 Results, Covariate selection: 'No other patient characteristics or laboratory values ... significantly affected any of the PK parameters (P > 0.05)'; Supplementary material lists age among the >10 covariates evaluated). Median (range) 14 months (3 months-13 years) overall; 7 (3-60) months without and 20 (4-156) months with haemofilter (Table 1)."
    ),
    HT = list(
      description = "Body height",
      units       = "cm",
      type        = "continuous",
      notes       = "Screened and not retained. Median (range) 65 (53-131) cm without and 82.5 (60-146.5) cm with haemofilter (Table 1)."
    ),
    BSA = list(
      description = "Body surface area",
      units       = "m^2",
      type        = "continuous",
      notes       = "Screened and not retained. Median (range) 0.35 (0.24-1.18) m^2 without and 0.45 (0.29-1.32) m^2 with haemofilter (Table 1). Distinct from the retained haemofilter membrane surface area (FILT_SA_MED / FILT_SA_LARGE), which is a device property, not a body size."
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate (bedside Schwartz equation)",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Screened and not retained on any parameter, including CLKRT (Butragueno-Laiseca 2025 Results, Covariate selection and Discussion: 'renal function expressed as eGFR was not found to have an impact on CLKRT, consistent with the results of Aulin et al.'). Median (range) 112 (29-283) mL/min in the non-CKRT group; not reported for CKRT patients (Table 1). The Limitations attribute the null result to the small sample size and to 27% of patients lacking cystatin C measurements."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "Screened and not retained. Median (range) 0.32 (0.1-2.46) mg/dL in the non-CKRT group; not reported for CKRT patients (Table 1)."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/dL",
      type        = "continuous",
      notes       = "Screened and not retained (Butragueno-Laiseca 2025 Discussion: 'albumin concentration did not show covariate effects in any of the model parameters despite the fact that teicoplanin is highly bound to that plasma protein'; the authors attribute this to the narrow albumin range, matching Aulin et al.). Median (range) 3.3 (2.1-4.2) g/dL without and 3.4 (2.5-4.2) g/dL with haemofilter (Table 1)."
    ),
    HGB = list(
      description = "Haemoglobin",
      units       = "g/dL",
      type        = "continuous",
      notes       = "Screened and not retained. Median (range) 9.7 (6.7-12.3) without and 9.8 (8.4-12.3) with haemofilter (Table 1). Distinct from the retained HCT, which enters the post-filter plasma-flow correction rather than a structural parameter."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 26L,
    n_studies      = 1L,
    n_observations = 173L,
    age_range      = "3 months to 13 years",
    age_median     = "14 months",
    weight_range   = "4.3-44 kg",
    weight_median  = "8 kg",
    race_ethnicity = "Not reported.",
    disease_state  = paste(
      "Critically ill children admitted to a paediatric intensive care unit and treated",
      "with teicoplanin, 12 of them undergoing continuous kidney replacement therapy in",
      "continuous venovenous haemodiafiltration modality. Most were in the post-operative",
      "period of a congenital cardiopathy (93.3% of the non-CKRT group and 91.7% of the",
      "CKRT group carried a cardiopathy diagnosis; 60% and 83.3% respectively were",
      "post-operative). Mechanical ventilation in 14/15 non-CKRT and 12/12 CKRT patients;",
      "PRISM III score median 6 (non-CKRT) and 8 (CKRT); mortality 20% and 41.7%.",
      "Indications were suspected catheter-related bacteraemia (40.7%),",
      "ventilator-associated pneumonia (25.9%), mediastinitis (18.5%), and sepsis or",
      "endocarditis risk in patients with extracorporeal devices (4 cases)",
      "(Supplementary material, Patient infections). No patient experienced",
      "nephrotoxicity or hepatotoxicity."
    ),
    dose_range     = paste(
      "Standard paediatric regimen: three loading doses of 10 mg/kg every 12 h in all",
      "patients, followed 24 h later by a maintenance dose of 10 mg/kg every 24 h in",
      "patients without a haemofilter or 3.3 mg/kg every 24 h in patients undergoing CKRT.",
      "Each dose given as an intravenous infusion over 5 min. Sampling generally began",
      "after at least three doses (after two doses in four patients), with blood drawn",
      "before and 1, 6, 12 and 24 h after the start of the infusion; in CKRT patients,",
      "samples were drawn simultaneously from the pre- and post-filter ports of the",
      "Prismaflex device. Median (range) 5.5 (4-13) samples per patient."
    ),
    regions        = "Hospital General Universitario Gregorio Maranon, Madrid, Spain (single centre).",
    renal_function = paste(
      "Non-CKRT group: eGFR (bedside Schwartz) median 112 mL/min, range 29-283 mL/min;",
      "serum creatinine median 0.32 mg/dL, range 0.1-2.46 mg/dL. CKRT group: eGFR and",
      "creatinine not reported (Table 1). CKRT prescription: blood flow median 50 mL/min",
      "(34-150), replacement fluid flow 95 mL/h (20-800), dialysate flow 300 mL/h",
      "(100-1800), extraction flow 45 mL/h (30-120), total effluent flow (CKRT dose)",
      "56 mL/kg/h (50-85); pre-filter citrate anticoagulation (18 mmol/L) in six patients.",
      "Filter surface areas: small 0.2 m^2 (n = 7), medium 0.6 m^2 (n = 3), large 1.2 m^2",
      "(n = 2)."
    ),
    notes          = paste(
      "Baseline demographics from Butragueno-Laiseca 2025 Table 1, reported separately for",
      "the 15 patients without and the 12 patients with a haemofilter; the group sizes sum",
      "to 27 rather than 26 because one patient started without CKRT and required CKRT",
      "after the fifth dose, contributing to both groups (Table 1 footnote a). Male sex 75%",
      "in the non-CKRT group and 40% in the CKRT group; no overall percentage is reported.",
      "Age, weight and height did not differ significantly between the groups (P > 0.1,",
      "Mann-Whitney U-test). 173 concentration measurements were modelled (72 plasma from",
      "non-CKRT patients, 51 pre-filter and 50 post-filter from CKRT patients); all were",
      "above the limit of quantification (1.0 mg/L) and ranged from 2.1 to 86.6 mg/L. Note",
      "that the Abstract states 172 measurements while the Results state 173, which is the",
      "figure consistent with the 72 + 51 + 50 breakdown; see the vignette Errata.",
      "Concentrations were assayed by validated LC-MS/MS (Supplementary material,",
      "Analytical method). NONMEM 7.4 with FOCE on log-transformed concentrations;",
      "parameter uncertainty by sampling importance resampling (PsN)."
    )
  )

  ini({
    # -----------------------------------------------------------------------
    # STRUCTURAL PARAMETERS -- Butragueno-Laiseca 2025 Table 2 (selected
    # model). The paper decomposes total elimination from the central
    # compartment into a renal arm (CLR) and a haemofilter arm (CLKRT); the
    # non-renal arm (CLNR) was not identifiable and was removed from the
    # model structure (Results, Base population model). The two retained arms
    # are mutually exclusive rather than additive: Table 2 footnote a states
    # that "CLR in patients with CKRT was estimated to be zero", and the
    # Supplementary material adds "thus, CL equals CLKRT" for CKRT patients.
    #
    #   lcl_renal         <- log(theta_CLR)   = log(0.169)   (Table 2)
    #   lcl_hemodialysis  <- log(theta_CLKRT) = log(0.119)   (Table 2)
    #
    # lcl_hemodialysis is the typical CKRT clearance for the SMALL (0.2 m^2)
    # haemofilter, which is the reference level of the theta_FILT covariate.
    # -----------------------------------------------------------------------
    lcl_renal        <- log(0.169); label("Typical renal clearance CLR at 8 kg, off CKRT (L/h)")             # Table 2: theta_CLR = 0.169 L/h (RSE 12%; SIR 95% CI 0.139-0.208)
    lcl_hemodialysis <- log(0.119); label("Typical haemofilter clearance CLKRT with a small 0.2 m2 filter (L/h)")  # Table 2: theta_CLKRT = 0.119 L/h (RSE 8%; SIR 95% CI 0.103-0.136)
    lvc              <- log(1.56);  label("Typical central volume of distribution V1 at 8 kg (L)")           # Table 2: theta_V1 = 1.56 L (RSE 17%; SIR 95% CI 1.03-1.98)
    lvp              <- log(3.03);  label("Typical peripheral volume of distribution V2 (L)")                # Table 2: theta_V2 = 3.03 L (RSE 30%; SIR 95% CI 2.10-4.53)
    lq               <- log(0.292); label("Typical inter-compartmental clearance CLD (L/h)")                 # Table 2: theta_CLD = 0.292 L/h (RSE 53%; SIR 95% CI 0.167-0.614)

    # -----------------------------------------------------------------------
    # ALLOMETRIC EXPONENTS -- Results, Covariate selection: "Body weight
    # impacted V1 and CLR (P < 0.01) significantly and was incorporated in the
    # model through the allometric expression using exponents of 1 and 0.75
    # for V1 and CLR, respectively." The exponents are structural values the
    # authors imposed, not estimates -- Table 2 reports no RSE or confidence
    # interval for either -- so both are encoded with fixed().
    #
    # The reference weight of 8 kg is printed in the Parameter model column of
    # Table 2 as (WGT/8) and equals the cohort median weight (Results, Patient
    # population). It is confirmed arithmetically by the Discussion: "The model
    # predicts an increase in 0.016 L/h and 0.19 L, respectively, for an
    # increment of 1 kg of body weight" -- d/dWGT of 0.169*(WGT/8)^0.75 at
    # WGT = 8 is 0.169*0.75/8 = 0.0158 L/h/kg, and d/dWGT of 1.56*(WGT/8) is
    # 1.56/8 = 0.195 L/kg.
    #
    # Body weight was NOT retained on V2 or CLD ("including body weight as a
    # covariate in V2 and CLD did not improve model fit").
    # -----------------------------------------------------------------------
    e_wt_cl_renal <- fixed(0.75); label("Allometric exponent on (WT/8) for CLR (unitless)")  # Results, Covariate selection: exponent 0.75 on CLR
    e_wt_vc       <- fixed(1);    label("Allometric exponent on (WT/8) for V1 (unitless)")   # Results, Covariate selection: exponent 1 on V1

    # -----------------------------------------------------------------------
    # HAEMOFILTER-SIZE EFFECT ON THE CKRT CLEARANCE ARM -- Table 2:
    # CLKRT = theta_CLKRT * theta_FILT, with theta_FILT taking one value per
    # filter surface area:
    #
    #   Small  (0.2 m^2) = 1     (reference)
    #   Medium (0.6 m^2) = 3.58  (RSE 13%; SIR 95% CI 2.78-4.40)
    #   Large  (1.2 m^2) = 5.04  (RSE 14%; SIR 95% CI 3.80-6.26)
    #
    # Encoded as a pair of mutually exclusive binary indicators applied as
    # e_filt_sa_med_cl_hemodialysis^FILT_SA_MED *
    # e_filt_sa_large_cl_hemodialysis^FILT_SA_LARGE so that the small filter
    # (both indicators 0) recovers the reference multiplier of 1. This is the
    # same power-encoding precedent used by Goti 2018 vancomycin and by
    # Patel 2011 fluconazole for binary multiplicative covariate factors.
    # -----------------------------------------------------------------------
    e_filt_sa_med_cl_hemodialysis   <- 3.58; label("Multiplicative factor on CLKRT for a medium 0.6 m2 filter (vs small)")  # Table 2: theta_FILT Med = 3.58
    e_filt_sa_large_cl_hemodialysis <- 5.04; label("Multiplicative factor on CLKRT for a large 1.2 m2 filter (vs small)")   # Table 2: theta_FILT Large = 5.04

    # -----------------------------------------------------------------------
    # INTER-INDIVIDUAL VARIABILITY -- Table 2, "IIV (%) (RSE, %)" column.
    #
    # Methods: "Inter-individual variability (IIV) was described with the
    # exponential model". Table 2 footnote: "IIV is expressed as coefficient
    # of variation (CV, %) calculated as sqrt(exp(omega^2) - 1) * 100, where
    # omega^2 corresponds to the variance of the random effects." Inverting
    # that relation gives omega^2 = log(1 + CV^2), the standard lognormal
    # conversion.
    #
    #   IIV CLR   = 30.4 % CV (RSE 32%), eta-shrinkage 33%
    #   IIV CLKRT = 14.1 % CV (RSE 19%), eta-shrinkage 45%
    #   IIV V1    = 39.2 % CV (RSE 45%), eta-shrinkage 40%
    #
    # No IIV on V2 or CLD: "Data supported the estimation of IIV for ... V1,
    # CLR, and CLKRT (P < 0.01), but not for V2 ... and CLD (P > 0.05)."
    # The etas are independent: "Covariance between random effects was also
    # non-significant (P > 0.05)."
    # -----------------------------------------------------------------------
    etalcl_renal        ~ log(1 + 0.304^2)  # Table 2: IIV CLR = 30.4 % CV
    etalcl_hemodialysis ~ log(1 + 0.141^2)  # Table 2: IIV CLKRT = 14.1 % CV
    etalvc              ~ log(1 + 0.392^2)  # Table 2: IIV V1 = 39.2 % CV

    # -----------------------------------------------------------------------
    # RESIDUAL ERROR -- Table 2, "Residual error, log (mg/mL)" rows.
    #
    # Methods: "For the analysis, teicoplanin concentrations were
    # logarithmically transformed. ... residual variability was described with
    # the additive error model in the logarithmic scale, and different
    # magnitudes of residual error were estimated for each type of measured
    # concentrations [plasma/pre-filter (CPre) and post-filter (CPost)]."
    #
    # An additive residual on log-transformed concentrations is a lognormal
    # residual on the linear concentration scale, so both magnitudes map to
    # the canonical expSd / lnorm() form rather than to propSd / addSd.
    #
    #   Plasma and pre-filter : 0.301 (RSE 8%;  SIR 95% CI 0.27-0.34), eps-shrinkage 9%
    #   Post-filter           : 0.333 (RSE 13%; SIR 95% CI 0.28-0.39)
    #
    # Plasma and pre-filter samples share one residual magnitude and are both
    # represented by the Cc observable; only the post-filter port has its own.
    # The "(mg/mL)" unit label in the Table 2 row header is a reporting slip --
    # a residual SD on the log scale is unitless, and every concentration in
    # the paper is in mg/L. See the vignette Errata.
    # -----------------------------------------------------------------------
    expSd              <- 0.301; label("Log-scale additive residual SD for plasma and pre-filter concentrations (unitless)")  # Table 2: plasma residual error 0.301
    expSd_Cpostfilter  <- 0.333; label("Log-scale additive residual SD for post-filter concentrations (unitless)")            # Table 2: post-filter residual error 0.333
  })

  model({
    # -------------------------------------------------------------------
    # 1. Individual PK parameters.
    #
    # The renal and haemofilter clearance arms are mutually exclusive, not
    # additive: RRT_CRRT_ACTIVE selects which one contributes. Off CKRT the
    # total clearance is the allometrically scaled renal clearance; on CKRT
    # the renal arm is zero (Table 2 footnote a) and the total clearance is
    # the filter-size-adjusted haemofilter clearance.
    # -------------------------------------------------------------------
    cl_renal <- (1 - RRT_CRRT_ACTIVE) * exp(lcl_renal + etalcl_renal) * (WT / 8)^e_wt_cl_renal

    cl_hemodialysis <-
      RRT_CRRT_ACTIVE * exp(lcl_hemodialysis + etalcl_hemodialysis) *
      e_filt_sa_med_cl_hemodialysis^FILT_SA_MED *
      e_filt_sa_large_cl_hemodialysis^FILT_SA_LARGE

    cl <- cl_renal + cl_hemodialysis
    vc <- exp(lvc + etalvc) * (WT / 8)^e_wt_vc
    vp <- exp(lvp)
    q  <- exp(lq)

    # 2. Micro-constants.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # 3. ODE system -- two-compartment IV disposition. Teicoplanin was given
    #    as a 5-min intravenous infusion (Methods, Dosing); the infusion
    #    duration is supplied on the dose records rather than modelled, since
    #    the paper does not estimate it.
    d/dt(central)     <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <-   k12 * central - k21 * peripheral1

    # -------------------------------------------------------------------
    # 4. Observations.
    #
    # Cc is the plasma concentration and, in CKRT patients, the pre-filter
    # concentration: the two are the same model prediction and share one
    # residual magnitude (Methods: residual error was estimated for
    # "plasma/pre-filter (CPre)" as one group and "post-filter (CPost)" as
    # the other).
    #
    # The post-filter concentration is the pre-filter concentration reduced
    # by the fraction of drug the haemofilter extracts in a single pass.
    # Supplementary material (Pharmacokinetic analysis):
    #
    #   Cpost = Cpre * (1 - CLKRT / phi_PlCorr)
    #
    # where phi_PlCorr is "the corrected plasma flow that goes into the
    # hemofilter, which is calculated as phi_Blood x BPR (blood to plasma
    # ratio)". The paper does not print the numerical definition of BPR;
    # it is encoded here as the plasma fraction of blood, 1 - HCT/100,
    # which is what makes the Results statement that "CKRT parameters and
    # haematocrit were included in the model" true of this equation. See
    # the vignette Errata.
    #
    # BFR is stored in the canonical mL/min and converted to L/h (x 0.06)
    # to match the L/h in which CLKRT is expressed.
    #
    # Both denominators carry a small numerical floor so that non-CKRT
    # subjects (BFR = 0, cl_hemodialysis = 0) and any extreme sampled
    # combination of a large filter with a low blood flow cannot produce a
    # division by zero or a negative concentration. Neither floor binds
    # anywhere in the covariate ranges the paper reports.
    # -------------------------------------------------------------------
    plasma_flow_filter <- max(BFR * 0.06 * (1 - HCT / 100), 0.001)
    filter_extraction  <- cl_hemodialysis / plasma_flow_filter

    Cc           <- central / vc
    Cpostfilter  <- Cc * max(1 - filter_extraction, 0.001)

    Cc          ~ lnorm(expSd)
    Cpostfilter ~ lnorm(expSd_Cpostfilter)
  })
}
