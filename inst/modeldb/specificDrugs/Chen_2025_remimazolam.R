Chen_2025_remimazolam <- function() {
  description <- paste(
    "Two-compartment population pharmacokinetic model for remimazolam",
    "given by continuous intravenous infusion for sedation in critically",
    "ill adults in the intensive care unit (Chen 2025). This is the first",
    "published popPK analysis of remimazolam in an ICU population, and",
    "the study's central finding is a negative one: of the 27 candidate",
    "covariates screened -- including body weight, age, sex, hepatic and",
    "renal function markers, and the presence of extracorporeal membrane",
    "oxygenation (ECMO) or continuous renal replacement therapy (CRRT)",
    "-- none survived backward elimination, so the final model carries no",
    "covariate effects at all. Disposition is described by first-order",
    "elimination from a central compartment plus a single peripheral",
    "compartment; inter-individual variability is exponential on all four",
    "structural parameters and residual error is purely proportional.",
    "The authors' companion Monte Carlo simulation shows that the",
    "context-sensitive half-time is essentially independent of infusion",
    "duration (15.6-21 min over infusions from 0.5 to 72 h), which is the",
    "property that makes the drug attractive for long-term ICU sedation.",
    "Note that the authors chose a two-compartment structure where most",
    "prior remimazolam analyses used three compartments; they attribute",
    "this to sparse arterial sampling confined to the post-infusion",
    "elimination phase (see the validation vignette).",
    sep = " "
  )
  reference <- paste(
    "Chen J, Wang X, Chen D, Liu X, Peng K, Tang R, Hu L, Wang Y, Bai Y,",
    "Chang L, Chen C. (2025). Population pharmacokinetic analysis of",
    "remimazolam after continuous infusion for sedation in critically ill",
    "patients. Front Pharmacol 16:1526266.",
    "doi:10.3389/fphar.2025.1526266.",
    sep = " "
  )
  vignette <- "Chen_2025_remimazolam"
  units    <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Remimazolam was given as an intravenous infusion by
  # micropump and quantified in plasma by HPLC-MS/MS (Chen 2025 Section 2.3),
  # so both states hold parent drug and the sampled matrix is plasma. The
  # inactive metabolite CNS 7054 was NOT modelled: Chen 2025 Section 4.4
  # states that the metabolite reference standard was unavailable because
  # remimazolam is a controlled drug, so no parent-metabolite model could be
  # built.
  compartmentData <- list(
    central     = list(analyte = "remimazolam", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "remimazolam", units = "mg", specimen = "plasma", verified = TRUE)
  )

  # No covariate was retained in the final model, so there is no
  # covariateData block. The 27 candidates below were each tested on CL and
  # on V1 by forward stepwise univariate inclusion (Chen 2025 Supplementary
  # Table S3, 54 numbered runs). Seven reached the forward-inclusion
  # threshold of dOFV <= -3.84, but none met the stricter backward
  # elimination criterion of dOFV >= +10.82, so the final model is the base
  # two-compartment model unchanged (Chen 2025 Section 3.2). No point
  # estimate is published for any of them, which is why they are documented
  # here rather than in covariateData.
  #
  # Six further screened covariates have no entry in the canonical covariate
  # register and are deliberately NOT given invented names here, since they
  # carry no point estimate and never appear in model(): ECMO status
  # (yes/no; runs 11-12), serum uric acid (runs 25-26), CKD-EPI estimated
  # glomerular filtration rate (runs 29-30), procalcitonin (runs 33-34),
  # arterial pH (runs 41-42) and platelet count (runs 43-44). They are listed
  # in full in the validation vignette's source-trace table.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age.",
      units       = "years",
      type        = "continuous",
      notes       = "Chen 2025 Supplementary Table S3 runs 5-6 (CL dOFV -5.564; V1 dOFV -4.610). Screened as a median-normalised exponential model (Chen 2025 Section 2.4, second displayed equation); not retained."
    ),
    SEXF = list(
      description = "Female sex indicator.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Chen 2025 Supplementary Table S3 runs 1-2 (CL dOFV -6.244, p < 0.01; V1 dOFV -6.037, p < 0.01). The strongest covariate signal in the screen, but it failed backward elimination. Chen 2025 Section 4.1 reports a post hoc simulation in which CL was 1.3-fold higher in females than in males (Supplementary Figure S3), and the authors discard the effect because men outnumbered women three to one (24 vs 8). Screened as a proportional model (Chen 2025 Section 2.4, first displayed equation) with the categorical covariate coded 0/1; the paper does not state which sex was coded 1, so the SEXF polarity here is the register default and is not asserted from the source."
    ),
    WT = list(
      description = "Total body weight at baseline.",
      units       = "kg",
      type        = "continuous",
      notes       = "Chen 2025 Supplementary Table S3 runs 3-4 (CL dOFV -4.610; V1 dOFV 0). Not retained. Chen 2025 Section 4.1 argues explicitly against weight-based dosing for remimazolam and notes that the 47-98 kg range was concentrated in 47-75 kg, with only three patients above it, so the negative result should not be extrapolated to obese patients."
    ),
    BMI = list(
      description = "Body mass index at baseline.",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Chen 2025 Supplementary Table S3 runs 7-8 (CL dOFV -4.610; V1 dOFV -4.610). Not retained."
    ),
    CREAT = list(
      description = "Serum creatinine.",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Chen 2025 Supplementary Table S3 runs 13-14 (CL dOFV -1.633; V1 dOFV -6.111, p < 0.01). Not retained. Reported in umol/L in Chen 2025 Table 1 (median 75.75, range 42.1-911.77)."
    ),
    CRCL = list(
      description = "Creatinine clearance by the Cockcroft-Gault equation.",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Chen 2025 Supplementary Table S3 runs 27-28 (CL dOFV 0; V1 dOFV 0) -- literally no improvement on either parameter. Not retained. Uncorrected mL/min (Cockcroft-Gault), NOT normalised to 1.73 m^2."
    ),
    BUN = list(
      description = "Blood urea nitrogen.",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Chen 2025 Supplementary Table S3 runs 53-54 (CL dOFV -3.678; V1 dOFV -4.610). Not retained. Chen 2025 Table 1 reports the unit as umol/L."
    ),
    CYSC = list(
      description = "Serum cystatin C.",
      units       = "mg/L",
      type        = "continuous",
      notes       = "Chen 2025 Supplementary Table S3 runs 49-50 (CL dOFV -4.810; V1 dOFV -4.610). Not retained."
    ),
    ALT = list(
      description = "Alanine aminotransferase.",
      units       = "U/L",
      type        = "continuous",
      notes       = "Chen 2025 Supplementary Table S3 runs 15-16. The V1 effect was the single largest in the whole screen (dOFV -11.367, p < 0.001), and Chen 2025 Section 4.1 reports that V1 rose with rising ALT, but it still failed backward elimination and was dropped. Supplementary Figure S3 shows the simulated concentration-time consequence."
    ),
    TPRO = list(
      description = "Total serum protein.",
      units       = "g/L",
      type        = "continuous",
      notes       = "Chen 2025 Supplementary Table S3 runs 17-18 (CL dOFV -6.441, p < 0.01; V1 dOFV 0). Not retained."
    ),
    ALB = list(
      description = "Serum albumin.",
      units       = "g/L",
      type        = "continuous",
      notes       = "Chen 2025 Supplementary Table S3 runs 19-20 (CL dOFV -4.610; V1 dOFV -5.110). Not retained."
    ),
    TBILI = list(
      description = "Total serum bilirubin.",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Chen 2025 Supplementary Table S3 runs 21-22 (CL dOFV -0.060; V1 dOFV 0). Not retained."
    ),
    DBIL = list(
      description = "Direct (conjugated) serum bilirubin.",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Chen 2025 Supplementary Table S3 runs 23-24 (CL dOFV -4.607; V1 dOFV 0). Not retained."
    ),
    CRP = list(
      description = "C-reactive protein.",
      units       = "mg/L",
      type        = "continuous",
      notes       = "Chen 2025 Supplementary Table S3 runs 51-52 (CL dOFV -6.752, p < 0.01; V1 dOFV 0). Not retained."
    ),
    WBC = list(
      description = "Peripheral leukocyte (white blood cell) count.",
      units       = "10^9/L",
      type        = "continuous",
      notes       = "Chen 2025 Supplementary Table S3 runs 31-32 (CL dOFV -4.610; V1 dOFV -4.610). Not retained."
    ),
    HGB = list(
      description = "Haemoglobin concentration.",
      units       = "g/L",
      type        = "continuous",
      notes       = "Chen 2025 Supplementary Table S3 runs 45-46 (CL dOFV -4.610; V1 dOFV -4.610). Not retained."
    ),
    HCT = list(
      description = "Haematocrit.",
      units       = "%",
      type        = "continuous",
      notes       = "Chen 2025 Supplementary Table S3 runs 47-48 (CL dOFV -4.610; V1 dOFV -4.610). Not retained. Relevant to ECMO circuit haemodilution, which the authors also screened."
    ),
    LACT = list(
      description = "Arterial lactate.",
      units       = "mmol/L",
      type        = "continuous",
      notes       = "Chen 2025 Supplementary Table S3 runs 35-36 (CL dOFV -4.610; V1 dOFV -7.671, p < 0.01). Not retained."
    ),
    SOD = list(
      description = "Arterial serum sodium.",
      units       = "mmol/L",
      type        = "continuous",
      notes       = "Chen 2025 Supplementary Table S3 runs 37-38 (CL dOFV -2.748; V1 dOFV -4.610). Not retained."
    ),
    POT = list(
      description = "Arterial serum potassium.",
      units       = "mmol/L",
      type        = "continuous",
      notes       = "Chen 2025 Supplementary Table S3 runs 39-40 (CL dOFV -4.610; V1 dOFV -5.314). Not retained."
    ),
    RRT_CRRT_STATUS = list(
      description = "Continuous renal replacement therapy in progress (1 = yes, 0 = no).",
      units       = "(binary)",
      type        = "binary",
      notes       = "Chen 2025 Supplementary Table S3 runs 9-10 (CL dOFV -0.226; V1 dOFV -0.120) -- among the weakest signals in the screen. 14 of 32 patients received CRRT (11 CVVH, 1 CVVHD, 2 CVVHDF). Chen 2025 Section 4.1 additionally reports that CRRT duration, blood flow rate and dialysate flow rate had no effect, and concludes that no dose adjustment is needed during CRRT."
    ),
    APACHE_II = list(
      description = "Acute Physiology and Chronic Health Evaluation II score on the day of ICU admission or drug administration.",
      units       = "points",
      type        = "continuous",
      notes       = "Recorded for the whole cohort (Chen 2025 Section 2.2; median 26, range 15-41 in Table 1) as the illness-severity descriptor, but it does not appear among the 54 numbered runs of Supplementary Table S3, so it was collected and reported rather than formally screened."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 32,
    n_studies      = 1,
    age_range      = "26-79 years",
    age_median     = "62 years",
    weight_range   = "47-98 kg",
    weight_median  = "63 kg",
    bmi_range      = "18.73-36.00 kg/m2",
    bmi_median     = "22.67 kg/m2",
    sex_female_pct = 25,
    race_ethnicity = c(Asian = 100),
    disease_state  = paste(
      "Critically ill adults requiring sedation in the intensive care unit.",
      "Admission diagnoses were respiratory failure (28.7%), heart failure",
      "(26.0%), severe pneumonia (15.1%), septic shock (12.3%), acute",
      "myocardial infarction (11.0%) and COPD (6.8%). Median APACHE II score",
      "26 (range 15-41); in-hospital mortality 25%. Organ function ranged",
      "from normal to severely impaired (ALT 7.8-549.2 U/L; creatinine",
      "clearance 8.69-193.83 mL/min). 18 of 32 patients were on ECMO",
      "(10 veno-arterial, 7 veno-venous, 1 veno-arterial-veno) and 14 of 32",
      "were on CRRT.",
      sep = " "
    ),
    dose_range     = paste(
      "Continuous intravenous infusion by micropump of remimazolam besylate",
      "at 1 or 2 mg/mL, titrated by the bedside clinician to the target",
      "Richmond Agitation-Sedation Scale score. Administered dose 2-17.28",
      "mg/h (median 6 mg/h); infusion duration 6.15-294.9 h (median 8.33 h).",
      "No loading dose was given as part of the observational protocol.",
      sep = " "
    ),
    regions        = "China (single centre; Maoming People's Hospital, Guangdong)",
    co_medication  = paste(
      "Remifentanil for analgesia in 96.8% of patients, which Chen 2025",
      "Section 4.3 invokes to explain why adequate sedation was achieved at",
      "plasma concentrations below the 400-1200 ng/mL range reported for",
      "remimazolam monotherapy. Also dexmedetomidine (12.5%), butorphanol",
      "(12.5%) and propofol (6.25%).",
      sep = " "
    ),
    n_observations = paste(
      "236 plasma concentrations from 243 collected (7 discarded as",
      "implausibly high or low). Arterial samples were drawn at the moment",
      "the infusion was stopped and at 10, 20, 30, 60, 90, 120 and 240 min",
      "afterwards, so the design informs the post-infusion elimination phase",
      "only. HPLC-MS/MS calibration range 1.0-1000 ng/mL.",
      sep = " "
    ),
    notes          = paste(
      "Prospective single-centre observational study run April-December",
      "2022 (ethics approval PJ2021MI-K009-01). Demographics reproduced from",
      "Chen 2025 Table 1. The cohort is described as Chinese by recruitment",
      "site; the paper does not tabulate race or ethnicity explicitly.",
      sep = " "
    )
  )

  ini({
    # -----------------------------------------------------------------
    # Structural parameters. Chen 2025 Table 2, "Estimate (% RSE)" column
    # of the two-compartment model, which is simultaneously the base and
    # the final model because no covariate was retained.
    #
    # Q is 20.0 L/h here, NOT the 21.9 L/h quoted in the Abstract and in
    # Section 3.2 -- 21.9 is the bootstrap MEDIAN (Table 2, "Median by
    # 1,000 bootstraps"), not the point estimate. Chen 2025 Supplementary
    # Table S4 settles it: the reported micro-constants k12 = 0.78 /h and
    # k21 = 0.58 /h require Q = 20.0 (20.0/25.5 = 0.784, 20.0/34.5 =
    # 0.580), whereas Q = 21.9 would give 0.859 and 0.635.
    # -----------------------------------------------------------------
    lcl <- log(58.2); label("Elimination clearance, CL (L/h)")                       # Chen 2025 Table 2: CL = 58.2 L/h (RSE 10.6%; bootstrap 95% CI 47.8-72.3)
    lvc <- log(25.5); label("Central compartment volume, V1 (L)")                    # Chen 2025 Table 2: V1 = 25.5 L (RSE 12.8%; bootstrap 95% CI 16.8-33.3)
    lq  <- log(20.0); label("Intercompartmental clearance, Q (L/h)")                 # Chen 2025 Table 2: Q = 20.0 L/h (RSE 21.8%; bootstrap median 21.9, 95% CI 12.2-34.55)
    lvp <- log(34.5); label("Peripheral compartment volume, V2 (L)")                 # Chen 2025 Table 2: V2 = 34.5 L (RSE 17.0%; bootstrap 95% CI 26.0-58.8)

    # -----------------------------------------------------------------
    # Inter-individual variability. Chen 2025 Section 2.4 specifies an
    # exponential IIV model (eta normally distributed, mean 0, variance
    # omega^2) and Table 2 reports each IIV as a percentage. The
    # percentages are read as 100 * omega, so the ini() variance is
    # (percentage/100)^2. That reading is what reconciles the two halves
    # of Table 2: the bootstrap columns are on the raw omega^2 scale, and
    # back-transforming them (omega = sqrt(median/100)) reproduces the
    # Estimate column to within 1.2% for the two well-identified IIVs --
    # CL (49.7 vs 50.3, shrinkage 1%) and Q (66.5 vs 65.7, shrinkage 27%).
    # The alternative reading omega^2 = log(CV^2 + 1) misses those same
    # two rows by 5% and 13.5%. See the vignette's "Assumptions and
    # deviations" section for the full round-trip table.
    # -----------------------------------------------------------------
    etalcl ~ 0.253009  # Chen 2025 Table 2: IIV-CL = 50.3% (RSE 15.6%, shrinkage 1%); 0.503^2
    etalvc ~ 0.027225  # Chen 2025 Table 2: IIV-V1 = 16.5% (RSE 68.7%, shrinkage 69%); 0.165^2
    etalq  ~ 0.431649  # Chen 2025 Table 2: IIV-Q  = 65.7% (RSE 50.0%, shrinkage 27%); 0.657^2
    etalvp ~ 0.374544  # Chen 2025 Table 2: IIV-V2 = 61.2% (RSE 24.7%, shrinkage 15%); 0.612^2

    # -----------------------------------------------------------------
    # Residual unexplained variability. Chen 2025 Section 2.4 says a
    # combined additive-plus-proportional model was used during model
    # building, but the final model retains the proportional term only:
    # Table 2 lists a single "Proportional error (%)" row with no additive
    # counterpart, and Section 3.2 states the two-compartment model was
    # selected "with a proportionality error of 25%".
    # -----------------------------------------------------------------
    propSd <- 0.25; label("Proportional residual error (fraction)")                  # Chen 2025 Table 2: proportional error = 25.0% (RSE 13.0%; bootstrap 95% CI 18.1-29.6)
  })

  model({
    # Individual parameters. No covariate enters the final model
    # (Chen 2025 Section 3.2), so these are the typical values times the
    # exponential IIV terms and nothing else.
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q  <- exp(lq  + etalq)
    vp <- exp(lvp + etalvp)

    # Micro-constants. Chen 2025 Supplementary Table S4 tabulates these
    # for the typical subject as k10 = 2.28 /h, k12 = 0.78 /h and
    # k21 = 0.58 /h, which the definitions below reproduce exactly.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment disposition with intravenous input into the central
    # compartment (Chen 2025 Section 2.2: micropump infusion; Section 3.2
    # selects the two-compartment structure).
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Amounts are carried in mg and volumes in L, so central/vc is mg/L.
    # Chen 2025 reports every plasma concentration in ng/mL (assay range
    # 1.0-1000 ng/mL, Section 2.3), so convert with 1 mg/L = 1000 ng/mL.
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
