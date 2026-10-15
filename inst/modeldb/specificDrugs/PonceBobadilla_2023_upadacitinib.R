PonceBobadilla_2023_upadacitinib <- function() {
  description <- paste0(
    "Two-compartment population PK model for oral upadacitinib (a selective ",
    "JAK1 inhibitor) re-estimated on a pooled dataset of phase I studies in ",
    "healthy volunteers, phase IIb studies in rheumatoid arthritis, Crohn's ",
    "disease, atopic dermatitis and ulcerative colitis, and the phase IIb/III ",
    "ulcerative colitis induction trials (1,927 subjects). Absorption is ",
    "formulation-dependent: the immediate-release capsule is absorbed ",
    "first-order from a depot after a lag time, while the extended-release ",
    "tablet uses a parallel mixed process in which 72.5% of the absorbed dose ",
    "enters the central compartment by a zero-order input of 3.39 h and the ",
    "remaining 27.5% is absorbed first-order from the depot after a lag time, ",
    "with 79.8% bioavailability relative to the immediate-release capsule. ",
    "Covariates: patient population (healthy volunteer and rheumatoid ",
    "arthritis versus atopic dermatitis / ulcerative colitis / Crohn's ",
    "disease), sex, creatinine clearance and aspartate aminotransferase on ",
    "CL/F, and sex and body weight on Vc/F. The proportional residual error ",
    "differs between phase I and phase II/III studies and all residual error ",
    "is inflated by a factor of 1.68 in the first two hours after a dose. ",
    "Successor to the cross-population model of Nader 2020; the companion ",
    "logistic exposure-response models of the same paper take the average ",
    "steady-state concentration from this model."
  )
  reference <- paste(
    "Ponce-Bobadilla AV, Stodtmann S, Eckert D, Zhou W, Liu W, Mohamed M-EF.",
    "Upadacitinib population pharmacokinetics and exposure-response",
    "relationships in ulcerative colitis patients. Clin Pharmacokinet.",
    "2023;62(1):101-112. doi:10.1007/s40262-022-01191-6.",
    "Parameter values are Supplementary Table 6 of the Online Supplementary",
    "Material (OSM); the covariate and residual-error equations are in the",
    "OSM Methods.",
    sep = " "
  )
  vignette <- "PonceBobadilla_2023_upadacitinib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Baseline total body weight.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Power-form effect on Vc/F only (exponent 0.493, OSM Table 6). The OSM ",
        "Methods state that 'continuous covariates were normalized to the ",
        "median of the population and incorporated in the model using a power ",
        "function' but the pooled-dataset median is not printed. A reference ",
        "of 70 kg is used: the rounded standard, matching the 70.35 kg median ",
        "body weight the paper prints for its ulcerative colitis induction ",
        "exposure-response population (OSM Table 8). See the vignette ",
        "Assumptions and deviations."
      ),
      source_name = "Body weight"
    ),
    SEXF = list(
      description = "Biological sex indicator (1 = female, 0 = male).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male). OSM Table 6 reports both sex effects as 'ratio of females compared to males'; the OSM Figure 1 forest plot uses male subjects as the reference group.",
      notes = paste0(
        "Multiplicative effects: CL/F x 0.874 and Vc/F x 0.811 in females ",
        "(OSM Table 6). The OSM Methods equation enters categorical covariates ",
        "as exp(cov * theta), which for a binary indicator is the tabulated ",
        "ratio raised to the indicator."
      ),
      source_name = "Sex"
    ),
    CRCL = list(
      description = paste0(
        "Baseline creatinine clearance in raw mL/min (NOT BSA-normalized), ",
        "as in the predecessor upadacitinib analyses."
      ),
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Power-form effect on CL/F (exponent 0.213, OSM Table 6). The pooled ",
        "median is not printed; the reference of 109 mL/min is the median the ",
        "paper prints for its maintenance-study PK population (OSM Table 4, ",
        "N = 623), and agrees with the 108.70 mL/min centering value of the ",
        "Klunder 2019 lineage (modellib('Klunder_2019_upadacitinib')). This is ",
        "the raw mL/min branch of the CRCL register entry, the same as the ",
        "sibling upadacitinib models: the OSM Figure 1 renal-function groups ",
        "(30-60, 60-90 and >= 90 mL/min) are in mL/min, and supplying a ",
        "BSA-normalized value would bias CL/F."
      ),
      source_name = "Creatinine clearance (CrCl)"
    ),
    AST = list(
      description = "Baseline serum aspartate aminotransferase.",
      units = "IU/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Power-form effect on CL/F (exponent -0.0810, OSM Table 6): higher AST ",
        "gives slightly lower clearance. The pooled median is not printed; a ",
        "reference of 18 IU/L is used, the rounded mean the paper prints for ",
        "its ulcerative colitis induction PK population (17.9 IU/L, OSM ",
        "Table 3). With an exponent of -0.081 a 20% error in the reference ",
        "moves CL/F by under 2%. The OSM Figure 1 forest plot compares AST >= ",
        "25 IU/L (the 90th percentile among ulcerative colitis patients) with ",
        "< 25 IU/L."
      ),
      source_name = "Aspartate amino transferase (AST)"
    ),
    DIS_HEALTHY = list(
      description = "Healthy-participant indicator (1 = healthy volunteer, 0 = patient).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (patient). OSM Table 6 reports the 'CL/F ratio of healthy patients compared to patients with AD, UC, or CD', so a patient with atopic dermatitis, ulcerative colitis or Crohn's disease is the reference.",
      notes = paste0(
        "Multiplicative effect on CL/F: healthy volunteers x 1.22 (OSM ",
        "Table 6). Set DIS_HEALTHY = 1 and DIS_RA = 0 for a healthy volunteer; ",
        "both are 0 for atopic dermatitis, ulcerative colitis and Crohn's ",
        "disease."
      ),
      source_name = "Subject population (HV vs. AD/UC/CD)"
    ),
    DIS_RA = list(
      description = "Rheumatoid arthritis indicator (1 = patient with rheumatoid arthritis, 0 = any other subject).",
      units = "(binary)",
      type = "binary",
      reference_category = "0. OSM Table 6 reports the 'CL/F ratio of RA patients compared to patients with AD, UC, or CD', so the reference is a patient with atopic dermatitis, ulcerative colitis or Crohn's disease.",
      notes = "Multiplicative effect on CL/F: rheumatoid arthritis x 0.929 (OSM Table 6).",
      source_name = "Subject population (RA vs. AD/UC/CD)"
    ),
    FORM_UPA_ER = list(
      description = "Formulation indicator for the dosing record (1 = extended-release tablet, 0 = immediate-release capsule).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (immediate-release capsule). OSM Table 6 reports the extended-release bioavailability relative to the immediate-release formulation (79.8%), so CL/F, Vc/F, Q/F and Vp/F are on immediate-release bioavailability.",
      notes = paste0(
        "Selects the absorption sub-model. Extended release: Ka 0.114 1/h, lag ",
        "0.159 h, 72.5% of the absorbed dose by a 3.39 h zero-order input and ",
        "79.8% relative bioavailability. Immediate release: Ka 1.85 1/h, lag ",
        "0.237 h, no zero-order arm. Each formulation carries its own Ka random ",
        "effect (OSM Table 6). All ulcerative colitis trials used the ",
        "extended-release tablet (OSM Table 1); the immediate-release capsule ",
        "appears only in the phase I and early phase II data."
      ),
      source_name = "Formulation (immediate release / extended release)"
    ),
    STUDY_UPA_PHASE1 = list(
      description = "Study-stratum indicator (1 = record from a phase I study, 0 = record from a phase II/III study).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (phase II/III studies, which supply all patient data including the ulcerative colitis trials).",
      notes = paste0(
        "Selects the proportional residual-error magnitude: SD 0.261 for ",
        "phase I and 0.456 for phase II/III (OSM Table 6). The OSM Methods ",
        "write the residual model as C = Chat * (1 + eps_Phase) + eps_ADD with ",
        "eps_Phase 'by clinical development Phase, i.e., Phase 1 versus Phase ",
        "2/3'. Unlike modellib('Klunder_2019_upadacitinib') the ISV magnitudes ",
        "are not split by phase in this model."
      ),
      source_name = "Clinical development phase (Phase 1 vs Phase 2/3)"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "upadacitinib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "upadacitinib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "upadacitinib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 1927L,
    age_range = "Not tabulated for the pooled dataset; ulcerative colitis induction patients mean 42.9 years (SD 14.4; OSM Table 3).",
    weight_range = "Not tabulated for the pooled dataset; ulcerative colitis induction patients mean 73.1 kg (SD 18.3; OSM Table 3).",
    sex_female_pct = 38,
    race_ethnicity = "Ulcerative colitis induction patients (N = 978): White 69%, Black or African American 3%, Asian 26%, Other/Multiple 2% (OSM Table 3).",
    disease_state = paste0(
      "Pooled healthy volunteers (phase I) and patients with rheumatoid ",
      "arthritis, Crohn's disease, atopic dermatitis (phase IIb) and moderately ",
      "to severely active ulcerative colitis (phase IIb/III). 978 ulcerative ",
      "colitis patients from induction Studies 1, 2 and 4 contributed to the ",
      "estimation; 623 patients from maintenance Study 3 were used for ",
      "post hoc validation."
    ),
    dose_range = "Ulcerative colitis: 7.5, 15, 30 or 45 mg once daily extended-release tablet for 8-week induction; 15 or 30 mg once daily for 52-week maintenance. Phase I and II data include immediate-release capsule doses.",
    regions = "Multinational; 22% of ulcerative colitis induction patients from the US and 23% from Asian regions (OSM Table 3).",
    sex_female_pct_note = "38% female among ulcerative colitis induction (N = 978) and maintenance (N = 623) PK patients (OSM Tables 3-4); the pooled-dataset percentage is not reported.",
    crcl_ulcerative_colitis = "Mean 112 mL/min (SD 32.9) in induction patients (OSM Table 3); median 109 mL/min (range 42.0-290) in maintenance patients (OSM Table 4).",
    notes = paste0(
      "NONMEM 7.4.4, FOCE with interaction. All parameters of the Nader 2020 ",
      "cross-population model were re-estimated after adding the ulcerative ",
      "colitis induction data up to Week 8; covariates were re-screened by ",
      "stepwise forward inclusion (alpha 0.01) and backward elimination ",
      "(alpha 0.001), which added AST on CL/F. Concentrations from the WuXi ",
      "laboratory in China were reduced by 11% for cross-laboratory ",
      "comparability. LLOQ 0.05 ng/mL; the first BLQ value after a dose was ",
      "kept at LLOQ/2 (OSM Methods)."
    )
  )

  ini({
    # All values are OSM Table 6 ('Summary Statistics of Pharmacokinetic
    # Parameter Estimates (Upadacitinib Population PK Final Model)'); each row
    # carries a 95% confidence interval, so every parameter was estimated and
    # none is wrapped in fixed().

    # ---------------- Disposition (immediate-release bioavailability) ------
    lcl <- log(37.0)
    label("Apparent oral clearance CL/F at the reference covariates (L/h)") # OSM Table 6 'CL/F (L/h)' = 37.0 (95% CI 35.8-38.2)
    lvc <- log(169)
    label("Apparent central volume Vc/F at the reference covariates (L)") # OSM Table 6 'Vc/F (L)' = 169 (95% CI 162-177)
    lq <- log(3.86)
    label("Apparent intercompartmental clearance Q/F (L/h)") # OSM Table 6 'Q/F (L/h)' = 3.86 (95% CI 3.69-4.03); the row label is merged with the tail of the preceding bioavailability row in the published table
    lvp <- log(61.6)
    label("Apparent peripheral volume Vp/F (L)") # OSM Table 6 'Vp/F (L)' = 61.6 (95% CI 60.4-62.7)

    # ---------------- Extended-release absorption ---------------------------
    lka_er <- log(0.114)
    label("First-order absorption rate constant of the extended-release tablet (1/h)") # OSM Table 6 'Extended-Release KA (1/h)' = 0.114 (95% CI 0.102-0.128)
    ltlag_er <- log(0.159)
    label("First-order absorption lag time of the extended-release tablet (h)") # OSM Table 6 'Extended-Release Absorption Lag Time (h)' = 0.159 (95% CI 0.150-0.169)
    logitffo <- log(0.275 / 0.725)
    label("Logit of the fraction of the extended-release dose taking the first-order depot route (unitless)") # OSM Table 6 'Fraction of Extended-Release Dose Absorbed through Zero-Order Process (%)' = 72.5 (95% CI 71.8-73.2); first-order share 1 - 0.725 = 0.275
    ld2 <- log(3.39)
    label("Duration of the zero-order extended-release input into the central compartment (h)") # OSM Table 6 'Zero-Order Absorption Duration (h)' = 3.39 (95% CI 3.32-3.45)
    lfrel_er <- log(0.798)
    label("Bioavailability of the extended-release tablet relative to the immediate-release capsule (fraction)") # OSM Table 6 'Bioavailability of the Extended-Release Formulation Relative to the Immediate-Release Formulation (%)' = 79.8 (95% CI 78.6-80.9)

    # ---------------- Immediate-release absorption --------------------------
    lka_ir <- log(1.85)
    label("First-order absorption rate constant of the immediate-release capsule (1/h)") # OSM Table 6 'Immediate-Release KA (1/h)' = 1.85 (95% CI 1.66-2.06)
    ltlag_ir <- log(0.237)
    label("Absorption lag time of the immediate-release capsule (h)") # OSM Table 6 'Immediate-Release Absorption Lag Time (h)' = 0.237 (95% CI 0.233-0.240)

    # ---------------- Covariate effects -------------------------------------
    e_crcl_cl <- 0.213
    label("Power exponent of creatinine clearance on CL/F (unitless)") # OSM Table 6 'Covariate Exponent of Creatinine Clearance on CL/F' = 0.213 (95% CI 0.158-0.268)
    e_ast_cl <- -0.0810
    label("Power exponent of aspartate aminotransferase on CL/F (unitless)") # OSM Table 6 'Covariate Exponent of Aspartate Amino Transferase on CL/F' = -0.0810 (95% CI -0.122 to -0.0402)
    e_sexf_cl <- log(0.874)
    label("Log ratio of CL/F in females relative to males (unitless)") # OSM Table 6 'CL/F Ratio of Females Compared to Males' = 0.874 (95% CI 0.837-0.913)
    e_dis_healthy_cl <- log(1.22)
    label("Log ratio of CL/F in healthy volunteers relative to AD/UC/CD patients (unitless)") # OSM Table 6 'CL/F Ratio of Healthy Patients Compared to Patients with AD, UC, or CD' = 1.22 (95% CI 1.13-1.32)
    e_dis_ra_cl <- log(0.929)
    label("Log ratio of CL/F in rheumatoid arthritis patients relative to AD/UC/CD patients (unitless)") # OSM Table 6 'CL/F Ratio of RA Patients Compared to Patients with AD, UC, or CD' = 0.929 (95% CI 0.888-0.972)
    e_wt_vc <- 0.493
    label("Power exponent of body weight on Vc/F (unitless)") # OSM Table 6 'Covariate Exponent of Body Weight on Vc/F' = 0.493 (95% CI 0.379-0.607)
    e_sexf_vc <- log(0.811)
    label("Log ratio of Vc/F in females relative to males (unitless)") # OSM Table 6 'Vc/F Ratio of Females Compared to Males' = 0.811 (95% CI 0.766-0.858)

    # ---------------- Inter-subject variability -----------------------------
    # OSM Table 6 gives each ISV as a percentage. Read here as
    # 100 * sqrt(omega^2), the convention the same group footnotes in the
    # predecessor upadacitinib analysis (Klunder 2019 Table 3: '%ISV was
    # calculated as SQRT(omega^2) x 100'), so omega^2 = (ISV/100)^2; the
    # printed 95% CIs, squared, are symmetric about omega^2 as a Wald interval
    # on the variance would be (immediate-release Ka: 0.629^2 and 0.790^2
    # average 0.5099 vs 0.714^2 = 0.5098). The OSM
    # Results state that the KA-CL/F and KA-Vc/F correlations were removed;
    # no CL/F-Vc/F covariance is reported, so the block is diagonal.
    etalcl ~ 0.119025 # OSM Table 6 'ISV on CL/F (%)' = 34.5; 0.345^2 = 0.119025
    etalvc ~ 0.0729 # OSM Table 6 'ISV on Vc/F (%)' = 27.0; 0.270^2 = 0.0729
    etalka_er ~ 0.889249 # OSM Table 6 'ISV on Extended-Release KA (%)' = 94.3; 0.943^2 = 0.889249
    etalka_ir ~ 0.509796 # OSM Table 6 'ISV on Immediate-Release KA (%)' = 71.4; 0.714^2 = 0.509796

    # ---------------- Residual error ----------------------------------------
    propSdPhase1 <- 0.261
    label("Proportional residual error SD, phase I studies (fraction)") # OSM Table 6 'Proportional Error (Phase 1) SD' = 0.261 (95% CI 0.258-0.264)
    propSdPhase23 <- 0.456
    label("Proportional residual error SD, phase II/III studies (fraction)") # OSM Table 6 'Proportional Error (Phase 2/3) SD' = 0.456 (95% CI 0.447-0.465)
    addSdBase <- 0.115
    label("Additive residual error SD (ng/mL)") # OSM Table 6 'Additive Error SD (ng/mL)' = 0.115 (95% CI 0.110-0.120)
    ruv_early <- 1.68
    label("Multiplicative scaling of the residual error within 2 h after a dose (unitless)") # OSM Table 6 'Scaling Factor on Residual Error' = 1.68 (95% CI 1.62-1.74); OSM Supplementary Results: 'a separate error term was introduced for the first two hours after dosing'
  })

  model({
    # Reference covariate values. The OSM Methods centre continuous
    # covariates on the pooled-population median, which is not printed; these
    # are the paper-printed anchors described in covariateData.
    ref_wt <- 70 # kg
    ref_crcl <- 109 # mL/min
    ref_ast <- 18 # IU/L

    # Formulation-gated absorption. For the immediate-release capsule there is
    # no zero-order arm, so ffo collapses to 1 and the whole dose goes through
    # the depot; exp(0) = 1 recovers unit relative bioavailability.
    ka_er <- exp(lka_er + etalka_er)
    ka_ir <- exp(lka_ir + etalka_ir)
    ka <- ka_er * FORM_UPA_ER + ka_ir * (1 - FORM_UPA_ER)
    tlag <- exp(ltlag_er * FORM_UPA_ER + ltlag_ir * (1 - FORM_UPA_ER))
    frel <- exp(lfrel_er * FORM_UPA_ER)
    ffo <- expit(logitffo) * FORM_UPA_ER + (1 - FORM_UPA_ER)

    # Individual parameters (OSM Methods: power functions of median-normalized
    # continuous covariates and exp(cov * theta) for binary covariates).
    cl <- exp(
      lcl + etalcl +
        e_crcl_cl * log(CRCL / ref_crcl) +
        e_ast_cl * log(AST / ref_ast) +
        e_sexf_cl * SEXF +
        e_dis_healthy_cl * DIS_HEALTHY +
        e_dis_ra_cl * DIS_RA
    )
    vc <- exp(lvc + etalvc + e_wt_vc * log(WT / ref_wt) + e_sexf_vc * SEXF)
    q <- exp(lq)
    vp <- exp(lvp)
    d2 <- exp(ld2)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Each EXTENDED-RELEASE administration is two dose records at the same
    # time and nominal amount: one to cmt = "depot" (first-order arm) and one
    # to cmt = "central" with rate = -2 (zero-order arm over the modelled
    # dur(central)). An IMMEDIATE-RELEASE administration needs only the depot
    # record because ffo = 1 makes f(central) zero.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - (kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(depot) <- frel * ffo
    f(central) <- frel * (1 - ffo)
    dur(central) <- d2
    # The paper describes the lag time as applying to the first-order
    # absorption, so the zero-order arm starts at the dose time.
    alag(depot) <- tlag

    # Doses in mg and volumes in L give mg/L; x 1000 converts to ng/mL.
    Cc <- 1000 * central / vc

    # Residual error: proportional SD by development phase plus a common
    # additive SD, both inflated by ruv_early within 2 h after the last dose.
    # tad(depot) is used because every dose, of either formulation, has a
    # depot record.
    early <- 0
    if (tad(depot) < 2) early <- 1
    ruv_scale <- 1 + (ruv_early - 1) * early
    propSd <- ruv_scale * (propSdPhase1 * STUDY_UPA_PHASE1 + propSdPhase23 * (1 - STUDY_UPA_PHASE1))
    addSd <- ruv_scale * addSdBase
    Cc ~ add(addSd) + prop(propSd)
  })
}
