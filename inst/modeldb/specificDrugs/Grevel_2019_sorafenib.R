Grevel_2019_sorafenib <- function() {
  description <- "One-compartment population PK model for oral sorafenib pooled from 10 phase I-IV studies in healthy volunteers and patients with hepatocellular, renal cell and differentiated thyroid carcinoma (Grevel 2019). Absorption is a chain of three sequential transit compartments sharing one first-order rate constant ktr (the dose enters the first; the third empties into the central compartment at ktr), and elimination is linear. Apparent clearance is 21.1% lower in women and 50.6% lower in the DECISION trial (radioiodine-refractory differentiated thyroid cancer), effects applied multiplicatively as (1 + theta * indicator). Relative bioavailability is 1 at steady state in the patient studies and 1.9 in the single-dose phase I healthy-volunteer studies. IIV on CL/F and V/F; proportional residual error whose magnitude differs between occasions with several samples (0.3) and occasions with a single sample (0.5). The paper's parametric time-to-event models for progression-free survival and adverse events are not included because their parameter estimates are not reported."
  reference <- paste(
    "Grevel J, Jentsch G, Austin R, Prins NH, Lettieri J, Mitchell D,",
    "Huang F, Brose MS, Schlumberger M, Meinhardt G, Pena CEA, Ploeger BA.",
    "Exposure-Response Modeling and Simulation of Progression-Free Survival",
    "and Adverse Events of Sorafenib Treatment in Patients With Advanced",
    "Thyroid Cancer.",
    "Clin Transl Sci. 2019;12(5):459-469.",
    "doi:10.1111/cts.12634.",
    sep = " "
  )
  vignette <- "Grevel_2019_sorafenib"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(analyte = "sorafenib", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "sorafenib", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "sorafenib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "sorafenib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    SEXF = list(
      description = "Sex indicator, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Figure S2: CL_cov,gender = (1 + theta8 * (Gender == female)); 21.1% lower CL/F in women (Table 1).",
      source_name = "Gender"
    ),
    STUDY_DECISION = list(
      description = "DECISION trial indicator (study 14295, NCT00984282; radioiodine-refractory differentiated thyroid cancer), 1 = DECISION, 0 = any other pooled study",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (the HCC / RCC patient studies and the healthy-volunteer studies)",
      notes = "Figure S2: CL_cov,study14295 = (1 + theta9 * (Study == 14295)); 50.6% lower CL/F (Table 1). The paper states the study effect cannot be separated from an effect of the DTC indication, because DECISION was the only DTC study.",
      source_name = "Study == 14295"
    ),
    STUDY_SORAFENIB_PHASE1 = list(
      description = "Single-dose phase I healthy-volunteer study indicator (studies 100483 and 100484), 1 = phase I single-dose study, 0 = patient study dosed to steady state",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (patient studies at steady state, F1 = 1)",
      notes = "Figure S2: F1 = 1 at steady state and theta7 for single-dose administration; Table S3 labels theta7 'F1, phase 1'. In this dataset the single-dose records are exactly the two phase I healthy-volunteer studies (Table S1).",
      source_name = "Phase 1 / single-dose flag"
    ),
    SAMPLE_INTENSIVE = list(
      description = "Per-observation indicator, 1 = the sample belongs to an occasion with several PK samples, 0 = an occasion with a single PK sample",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (single-observation occasion)",
      notes = "Figure S2: sigma_nSamples = theta4 for multiple-observation occasions and theta6 for single-observation occasions. Single-sample occasions are the sparse patient samples (e.g. the one unsupervised sample per DECISION patient); multiple-sample occasions include the dense phase I profiles.",
      source_name = "nSamples"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened (Methods, PK modeling) but not retained in the final model."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened but not retained."
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened (Asian vs non-Asian); no significant effect (Results)."
    ),
    AGE = list(description = "Age", units = "years", type = "continuous", notes = "Screened but not retained."),
    AST = list(
      description = "Aspartate aminotransferase (SGOT)",
      units = "U/L",
      type = "continuous",
      notes = "Screened but not retained; no significant hepatic-function effect (Results)."
    ),
    ALT = list(
      description = "Alanine aminotransferase (SGPT)",
      units = "U/L",
      type = "continuous",
      notes = "Screened but not retained."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "mg/dL",
      type = "continuous",
      notes = "Screened but not retained."
    ),
    ALP = list(
      description = "Alkaline phosphatase",
      units = "U/L",
      type = "continuous",
      notes = "Screened but not retained."
    ),
    LDH = list(
      description = "Lactate dehydrogenase",
      units = "U/L",
      type = "continuous",
      notes = "Screened but not retained."
    ),
    INR_BASE = list(
      description = "Prothrombin time international normalized ratio",
      units = "(ratio)",
      type = "continuous",
      notes = "Screened but not retained."
    ),
    CRCL = list(
      description = "Creatinine clearance",
      units = "mL/min",
      type = "continuous",
      notes = "Screened but not retained; no significant renal-function effect (Results)."
    ),
    CONMED_CYP3A4_IND = list(
      description = "Concomitant CYP3A4 inducer",
      units = "(binary)",
      type = "binary",
      notes = "Screened but not retained."
    ),
    CONMED_CYP3A4_INH = list(
      description = "Concomitant CYP3A4 inhibitor",
      units = "(binary)",
      type = "binary",
      notes = "Screened but not retained."
    ),
    CONMED_UGT1A9_INH = list(
      description = "Concomitant UGT1A9 inhibitor",
      units = "(binary)",
      type = "binary",
      notes = "Screened but not retained."
    ),
    CONMED_UGT1A9_IND = list(
      description = "Concomitant UGT1A9 inducer",
      units = "(binary)",
      type = "binary",
      notes = "Screened but not retained."
    ),
    CONMED_LEVOTHYROXINE = list(
      description = "Concomitant levothyroxine",
      units = "(binary)",
      type = "binary",
      notes = "Screened but not retained."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 859L,
    n_studies = 10L,
    n_observations = 3141L,
    age_range = "21-97 years",
    age_median = "62 years",
    weight_range = "34-140 kg",
    weight_median = "70 kg",
    sex_female_pct = 24.6,
    race_ethnicity = c(Asian = 35, NonAsian = 65),
    disease_state = "Hepatocellular carcinoma (n = 332), renal cell carcinoma (n = 332), radioiodine-refractory differentiated thyroid cancer (DECISION, n = 156) and healthy volunteers (n = 39).",
    dose_range = "Single oral doses in the phase I healthy-volunteer studies; 400 mg twice daily continuous (with protocol dose reductions to 600, 400 and 200 mg/day) in the patient studies.",
    regions = "Multinational (SHARP, TARGET, DECISION; Asia-Pacific and Japanese phase II/III studies; Taiwanese phase IV HATT study).",
    notes = "Demographics from Table S2 (median (range)); studies from Table S1. Dense PK in the 39 healthy volunteers, sparse PK in patients; one unsupervised PK sample per DECISION patient. Sorafenib concentrations below the LLOQ (1-10 ug/L by assay) were excluded."
  )

  ini({
    lktr <- log(1.0)
    label("Transit rate constant ktr (1/h)") # Table S3 theta3 ktr = 1.0 1/h (95% CI 0.9-1.2)
    lcl <- log(9.4)
    label("Apparent clearance CL/F, male, non-DECISION study (L/h)") # Table S3 theta1 CL/F = 9.4 L/h (95% CI 8.9-9.9); Table 1 male HCC/RCC CL/F 9.4 L/h
    lvc <- log(554)
    label("Apparent volume of distribution V1/F (L)") # Table S3 theta2 V1/F = 554.0 L (95% CI 440-668)
    lfdepot <- fixed(log(1))
    label("Relative bioavailability at steady state in patient studies (unitless)") # Figure S2 'F1 = 1 steady state'

    e_study_phase1_fdepot <- 1.9
    label("Relative bioavailability F1 in the single-dose phase I studies (unitless)") # Table S3 theta7 'F1, phase 1' = 1.9 (95% CI 1.6-2.3); Figure S2 'F1 = theta7 single-dose'
    e_sexf_cl <- -0.211
    label("Fractional change in CL/F for women (unitless)") # Table 1 'Difference, female vs. male' CL/F = -21.1%; Table S3 theta8 = -0.2 (95% CI -0.3 to -0.1)
    e_study_decision_cl <- -0.506
    label("Fractional change in CL/F in the DECISION trial (unitless)") # Table 1 'Difference, DTC vs. HCC, RCC' CL/F = -50.6%; Table S3 theta9 = -0.5 (95% CI -0.6 to -0.5)

    etalcl ~ 0.2 # Table S3 omega(1,1) IIV CL = 0.2 (95% CI 0.2-0.3)
    etalvc ~ 0.4 # Table S3 omega(2,2) IIV V1 = 0.4 (95% CI 0.2-0.7)

    propSdIntensive <- 0.3
    label("Proportional residual SD, occasions with several samples (fraction)") # Table S3 theta4 'Proportional error' = 0.3 (95% CI 0.3-0.4); sigma(1,1) = 1 fixed
    propSdSparse <- 0.5
    label("Proportional residual SD, occasions with a single sample (fraction)") # Table S3 theta6 'Proportional error single-observation' = 0.5 (95% CI 0.5-0.6)
  })

  model({
    # Figure S2: CL_i = CL_TV * CL_cov * exp(eta1),
    # CL_cov = (1 + theta8 * female) * (1 + theta9 * (Study == 14295))
    ktr <- exp(lktr)
    cl <- exp(lcl + etalcl) *
      (1 + e_sexf_cl * SEXF) *
      (1 + e_study_decision_cl * STUDY_DECISION)
    vc <- exp(lvc + etalvc)
    kel <- cl / vc

    # Figure S2: A1-A3 are the three transit compartments sharing ktr; the
    # dose enters A1 (depot) and A3 (transit2) empties into A4 (central).
    d / dt(depot) <- -ktr * depot
    d / dt(transit1) <- ktr * depot - ktr * transit1
    d / dt(transit2) <- ktr * transit1 - ktr * transit2
    d / dt(central) <- ktr * transit2 - kel * central

    # Figure S2: F1 = 1 at steady state (patient studies), theta7 for the
    # single-dose phase I studies.
    fdepot <- exp(lfdepot) *
      (1 + (e_study_phase1_fdepot - 1) * STUDY_SORAFENIB_PHASE1)
    f(depot) <- fdepot

    Cc <- central / vc

    # Figure S2: Ct = A4 / V1 * (1 + eps), eps ~ N(0, sigma_nSamples^2) with
    # sigma = theta4 on multiple-observation and theta6 on single-observation
    # occasions; the additive component (theta5) is fixed to 0 and omitted.
    propSd <- propSdIntensive * SAMPLE_INTENSIVE +
      propSdSparse * (1 - SAMPLE_INTENSIVE)
    Cc ~ prop(propSd)
  })
}
