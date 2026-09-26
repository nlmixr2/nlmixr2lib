Lin_2020_glasdegib <- function() {
  description <- paste(
    "Two-compartment population PK model with first-order absorption for",
    "oral glasdegib (a Hedgehog-pathway SMO inhibitor) in 269 adults with",
    "advanced hematologic malignancies (AML, MDS and others; studies",
    "B1371001 and B1371003) or solid tumors (study B1371002) given 5-640 mg",
    "once daily (Lin 2020). Allometric baseline body weight scaling",
    "(exponent 0.75 on CL/F and Q/F, 1 on Vc/F and Vp/F, fixed, 70 kg",
    "reference). CL/F additionally depends on weight-standardized",
    "Cockcroft-Gault creatinine clearance (power), baseline percentage bone",
    "marrow blasts (linear, centered at 38.2%) and concomitant moderate or",
    "strong CYP3A inhibitors (linear); a solid-tumor indicator lowers Vp/F",
    "and Q/F. Exponential IIV on CL/F, Vc/F, Vp/F, Q/F and ka. Residual",
    "error is additive on the log scale with separate SDs for hematologic",
    "and solid-tumor patients.",
    sep = " "
  )
  reference <- paste(
    "Lin S, Shaik N, Martinelli G, Wagner AJ, Cortes J, Ruiz-Garcia A.",
    "Population Pharmacokinetics of Glasdegib in Patients With Advanced",
    "Hematologic Malignancies and Solid Tumors.",
    "J Clin Pharmacol. 2020;60(5):605-616. doi:10.1002/jcph.1556.",
    "Studies B1371001 (NCT00953758), B1371002 (NCT01286467) and",
    "B1371003 (NCT01546038).",
    sep = " "
  )
  vignette <- "Lin_2020_glasdegib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Two residual SDs selected by tumor type (Results 'Base Model': '2
  # thetarized residual proportional errors separately for hematology and
  # solid tumor patients').
  paper_specific_residual_sds <- c("expSdHeme", "expSdSolid")

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Baseline (time-fixed) weight. Allometric scaling referenced to",
        "70 kg with fixed exponents 0.75 (CL/F, Q/F) and 1 (Vc/F, Vp/F)",
        "(Results 'Base Model' and the final-model equations). Also used",
        "to weight-standardize creatinine clearance (see CRCL)."
      ),
      source_name = "BWT"
    ),
    CRCL = list(
      description = "Baseline creatinine clearance by the Cockcroft-Gault equation, NOT size-normalized",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Supply the RAW Cockcroft-Gault value in mL/min (neither",
        "BSA-normalized nor per-70-kg). The paper's covariate is the",
        "weight-standardized value WNCL = BCCL x 70 / BWT (Table 1",
        "footnote a), which model() computes from CRCL and WT; it enters",
        "CL/F as (WNCL / 71.23)^0.406, 71.23 mL/min being the median",
        "WNCL. The paper standardized to 70 kg to avoid counting the",
        "body-weight effect twice alongside the allometric term. Cohort",
        "median raw CRCL 80.9 mL/min, range 31.4-238.4 (Table 2)."
      ),
      source_name = "BCCL"
    ),
    BMBLAST_PCT = list(
      description = "Baseline percentage of blasts in bone marrow",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Linear effect on CL/F centered at 38.2% (the median among the",
        "hematologic-malignancy patients). Not collected in solid-tumor",
        "patients (Table 2: 'Not applicable'); the paper imputes a",
        "missing continuous covariate at the population median (Methods",
        "'Covariate Analyses'), so set BMBLAST_PCT = 38.2 for solid-tumor",
        "patients, which makes the term equal 1. Hematologic cohort",
        "median 39.3%, range 0-100 (Table 2)."
      ),
      source_name = "BPBL"
    ),
    CONMED_CYP3A4_INH_MOD = list(
      description = "Concomitant moderate CYP3A inhibitor, 1 = yes, 0 = no",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant moderate CYP3A inhibitor)",
      notes = paste(
        "Linear effect on CL/F, (1 - 0.173) = 17% lower CL/F. The paper",
        "does not say whether the flag was per record or per subject;",
        "Table 2 summarises one (most extreme) record per patient: 62 of",
        "272 patients moderate. The bootstrap 95% CI of the coefficient",
        "includes 0 (Table 3). A subject with both flags set gets the",
        "product of the two terms, as written in the paper's equation."
      ),
      source_name = "CYPmoderate"
    ),
    CONMED_CYP3A4_INH_STRONG = list(
      description = "Concomitant strong CYP3A inhibitor, 1 = yes, 0 = no",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant strong CYP3A inhibitor)",
      notes = paste(
        "Linear effect on CL/F, (1 - 0.303) = 30% lower CL/F. 56 of 272",
        "patients (Table 2). Mostly azole antifungals in the AML / MDS",
        "patients of B1371003 (Discussion)."
      ),
      source_name = "CYPstrong"
    ),
    TUMTP_SOLID = list(
      description = "Solid-tumor malignancy, 1 = solid tumor, 0 = hematologic malignancy",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (hematologic malignancy)",
      notes = paste(
        "Linear effects on Vp/F (1 - 0.825) and Q/F (1 - 0.653); also",
        "selects the solid-tumor residual SD. The 23 patients of study",
        "B1371002 are the solid-tumor patients."
      ),
      source_name = "Solid"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Baseline age",
      units = "years",
      type = "continuous",
      notes = "Tested on CL/F (Table 1); not retained. Cohort median 69 years, range 25-92 (Table 2)."
    ),
    SEXF = list(
      description = "Female sex, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      notes = "Tested on CL/F and Vc/F (Table 1); not retained. 91 of 272 patients female (Table 2)."
    ),
    AST = list(
      description = "Baseline aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Tested on CL/F as BAST (Table 1); not retained. Median 21.0 U/L (Table 2)."
    ),
    TBILI = list(
      description = "Baseline total bilirubin",
      units = "mg/dL",
      type = "continuous",
      notes = "Tested on CL/F as BBIL (Table 1); not retained. Median 0.6 mg/dL (Table 2)."
    ),
    ALB = list(
      description = "Baseline serum albumin",
      units = "g/dL",
      type = "continuous",
      notes = "Tested on CL/F, Vc/F and Vp/F as BALB (Table 1); not retained. Median 3.7 g/dL (Table 2)."
    )
  )

  compartmentData <- list(
    depot = list(analyte = "glasdegib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "glasdegib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "glasdegib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 269L,
    n_studies = 3L,
    n_observations = "3616 plasma glasdegib concentrations; no sample was below the 0.200 ng/mL LLOQ",
    age_range = "25-92 years (median 69) (Table 2, all 272 enrolled patients)",
    weight_range = "43.5-145.6 kg (median 78.6) (Table 2)",
    sex_female_pct = 33.5,
    race_ethnicity = c(White = 86.4, Black = 5.9, Asian = 3.3, Hispanic = 4.0, Other = 0.4),
    disease_state = paste(
      "Advanced hematologic malignancies (246 patients; AML or high-risk",
      "MDS in B1371003, mixed hematologic malignancies in B1371001) and",
      "solid tumors (23 patients, B1371002)"
    ),
    dose_range = paste(
      "Oral glasdegib 5-640 mg once daily; 187 of 272 patients (69%) at",
      "the clinical dose of 100 mg QD (Supplemental Table S3). Fasted in",
      "B1371001 / B1371002, without regard to food in B1371003, where",
      "glasdegib was combined with low-dose cytarabine, decitabine or",
      "7 + 3 chemotherapy."
    ),
    renal_function = "CRCL median 80.9 mL/min (31.4-238.4); 106 normal, 103 mild, 62 moderate impairment, none severe (Table 2)",
    hepatic_function = "220 normal, 43 mild, 3 moderate, 1 severe impairment (NCI ODWG) (Table 2)",
    co_medication = "Moderate CYP3A inhibitor 62 patients, strong CYP3A inhibitor 56 patients (Table 2)",
    notes = paste(
      "Demographics in Table 2 are for the 272 enrolled patients; 3",
      "without PK samples were excluded, leaving 269 in the analysis.",
      "NONMEM 7.3, SAEM followed by importance sampling; 1000-replicate",
      "bootstrap."
    )
  )

  ini({
    # Final-model typical values: Lin 2020 Table 3 'Final Model Estimate'
    # column; covariate coefficients from the final-model equations in
    # Results 'Covariate Analyses and Final Model', which print more digits
    # than Table 3 (e.g. 0.173 vs -0.17).
    lka <- log(0.06); label("First-order absorption rate constant ka (1/h)") # Table 3 'ka, hour-1' = 0.06 (RSE 5.0%)
    lcl <- log(6.27); label("Apparent clearance CL/F for a 70 kg hematologic patient, WNCL 71.23 mL/min, 38.2% marrow blasts, no CYP3A inhibitor (L/h)") # Table 3 'CL/F, L/h' = 6.27 (RSE 6.4%); CL/F equation
    lvc <- log(3.32); label("Apparent central volume Vc/F at 70 kg (L)") # Table 3 'Vc/F, L' = 3.32 (RSE 23.7%); Vc/F equation
    lvp <- log(279.21); label("Apparent peripheral volume Vp/F at 70 kg, hematologic patient (L)") # Table 3 'Vp/F, L' = 279.21 (RSE 90.0%); Vp/F equation
    lq <- log(1.29); label("Apparent intercompartmental clearance Q/F at 70 kg, hematologic patient (L/h)") # Table 3 'Q/F, L/h' = 1.29 (RSE 46.6%); Q/F equation (Results prose: 1.288)

    e_wt_cl_q <- fixed(0.75); label("Allometric exponent of baseline weight on CL/F and Q/F (unitless)") # Results 'Base Model': scaling factor of 0.75 on CL/F and Q/F; Results 'Body Weight': 'fixed effect'
    e_wt_vc_vp <- fixed(1); label("Allometric exponent of baseline weight on Vc/F and Vp/F (unitless)") # Results 'Base Model': 1.0 on Vc/F and Vp/F

    e_crcl_cl <- 0.406; label("Power exponent of weight-standardized CRCL (WNCL/71.23) on CL/F (unitless)") # CL/F equation exponent 0.406; Table 3 'theta WNCL' = 0.41 (RSE 23.3%)
    e_bmblast_pct_cl <- -0.004; label("Linear effect of baseline % bone marrow blasts (centered at 38.2%) on CL/F (1/%)") # CL/F equation '(1 - 0.004 (BPBL - 38.20))'; Table 3 'theta BPBL' = -0.004 (RSE 34.5%)
    e_cyp3a4_inh_mod_cl <- -0.173; label("Fractional change in CL/F with a moderate CYP3A inhibitor (unitless)") # CL/F equation '(1 - 0.173 CYPmoderate)'; Table 3 'theta CYP mod' = -0.17 (RSE 59.4%)
    e_cyp3a4_inh_strong_cl <- -0.303; label("Fractional change in CL/F with a strong CYP3A inhibitor (unitless)") # CL/F equation '(1 - 0.303 CYPstrong)'; Table 3 'theta CYP strong' = -0.30 (RSE 36.6%)
    e_tumtp_solid_vp <- -0.825; label("Fractional change in Vp/F for solid-tumor patients (unitless)") # Vp/F equation '(1 - 0.825 Solid)'; Table 3 'theta solid' (Vp/F) = -0.825 (RSE 14.6%)
    e_tumtp_solid_q <- -0.653; label("Fractional change in Q/F for solid-tumor patients (unitless)") # Q/F equation '(1 - 0.653 Solid)'; Table 3 'theta solid' (Q/F) = -0.65 (RSE 24.4%)

    # IIV: Table 3 reports CV%; converted with omega^2 = log(1 + CV^2).
    etalcl ~ 0.169658 # Table 3 IIV 'CL/F' = 43.0% (shrinkage 16.2%)
    etalvc ~ 1.72907 # Table 3 IIV 'Vc/F' = 215.3% (shrinkage 22.6%)
    etalvp ~ 0.82083 # Table 3 IIV 'Vp/F' = 112.8% (shrinkage 69.1%)
    etalq ~ 0.359745 # Table 3 IIV 'Q/F' = 65.8% (shrinkage 62.9%)
    etalka ~ 0.0180609 # Table 3 IIV 'ka' = 13.5% (shrinkage 70.7%)

    # Residual error: log-transformed concentrations with a thetarized
    # sigma (Methods 'Structural Model Development'), i.e. additive on the
    # log scale with SD = theta.
    expSdHeme <- 0.658; label("Log-scale residual SD, hematologic-malignancy patients (unitless)") # Table 3 'Residual error, proportional error: Hematologic' = 65.8% (RSE 5.2%)
    expSdSolid <- 0.595; label("Log-scale residual SD, solid-tumor patients (unitless)") # Table 3 'Residual error, proportional error: Solid tumors' = 59.5% (RSE 7.0%)
  })
  model({
    # Weight-standardized creatinine clearance (Table 1 footnote a)
    wncl <- CRCL * 70 / WT

    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl_q *
      (1 + e_cyp3a4_inh_mod_cl * CONMED_CYP3A4_INH_MOD) *
      (1 + e_cyp3a4_inh_strong_cl * CONMED_CYP3A4_INH_STRONG) *
      (wncl / 71.23)^e_crcl_cl *
      (1 + e_bmblast_pct_cl * (BMBLAST_PCT - 38.2))
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc_vp
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vc_vp * (1 + e_tumtp_solid_vp * TUMTP_SOLID)
    q <- exp(lq + etalq) * (WT / 70)^e_wt_cl_q * (1 + e_tumtp_solid_q * TUMTP_SOLID)
    ka <- exp(lka + etalka)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # mg / L = ug/mL; x 1000 gives ng/mL
    Cc <- 1000 * central / vc
    expSdi <- expSdHeme * (1 - TUMTP_SOLID) + expSdSolid * TUMTP_SOLID
    Cc ~ lnorm(expSdi)
  })
}
