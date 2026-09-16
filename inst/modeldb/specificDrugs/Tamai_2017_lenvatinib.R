Tamai_2017_lenvatinib <- function() {
  description <- paste0(
    "Three-compartment population PK model for the multikinase inhibitor ",
    "lenvatinib, re-estimated on a 452-subject / 8761-observation pool of ",
    "8 phase 1 clinical pharmacology studies in healthy adults, 4 phase 1 ",
    "dose-finding studies in mixed solid tumors, and the phase 1/2 study ",
    "202 in advanced hepatocellular carcinoma (HCC) Child-Pugh class A ",
    "(Tamai 2017, NCT00946153). Sequential zero-order release into the ",
    "depot over a duration D1 followed by first-order absorption (Ka), ",
    "with linear elimination from the central compartment. Covariate ",
    "effects on CL/F: body weight (power 0.708, shared with Q1/F and ",
    "Q2/F), CYP3A4 inducers (+30 percent), CYP3A inhibitors ",
    "(-7.8 percent), healthy-subject cohort (+19 percent vs cancer ",
    "patients), and alkaline phosphatase above the upper limit of normal ",
    "(-14.8 percent). Body weight also acts on V1/F, V2/F and V3/F with ",
    "a separately estimated power of 1.08, and capsule relative ",
    "bioavailability is 0.867 versus the tablet reference. This is the ",
    "HCC-focused successor to Gupta_2016_lenvatinib.R: the same structural ",
    "model re-fitted after adding study 202, which is why the estimates ",
    "differ modestly and the albumin term of the Gupta model is absent. ",
    "The companion Tamai_2017_lenvatinib_teae_dosemod.R carries the ",
    "paper's exposure-response logistic model for early dose modification."
  )
  reference <- paste(
    "Tamai T, Hayato S, Hojo S, Suzuki T, Okusaka T, Ikeda K, Kumada H.",
    "Dose finding of lenvatinib in subjects with advanced hepatocellular",
    "carcinoma based on population pharmacokinetic and exposure-response",
    "analyses.",
    "J Clin Pharmacol. 2017;57(9):1138-1147.",
    "doi:10.1002/jcph.917.",
    sep = " "
  )
  vignette <- "Tamai_2017_lenvatinib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight at baseline.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Centred at 75 kg, the value written into every Tamai 2017 Table 2",
        "covariate equation as '(WGT/75)'; the Table 1 pooled median is",
        "75.1 kg (mean 75.8, SD 17.7, range 42.7-147.0), so 75 is the",
        "rounded median the authors adopted. Two SEPARATE powers are",
        "estimated, not one allometric constant shared across the model:",
        "theta_WGT1 = 0.708 on CL/F, Q1/F and Q2/F, and theta_WGT2 = 1.08",
        "on V1/F, V2/F and V3/F. Both are ESTIMATED here (RSE 6.58 and 5.42",
        "percent, bootstrap 95 percent intervals 0.538-0.886 and",
        "0.876-1.28), unlike Gupta_2016_lenvatinib.R where the same two",
        "exponents were held fixed at the canonical 0.75 and 1. Body weight",
        "is the covariate the paper's whole dose-finding argument rests on:",
        "because CL/F rises less than proportionally with weight while dose",
        "is flat, AUC at steady state falls as weight rises, which is what",
        "motivates the 8 mg (< 60 kg) versus 12 mg (>= 60 kg) starting-dose",
        "split. Only 10.5 percent of CL/F interindividual variability is",
        "explained by body weight (Tamai 2017 Discussion)."
      ),
      source_name = "WGT"
    ),
    ALP = list(
      description = "Serum alkaline phosphatase activity.",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters the model ONLY as the binary indicator 'ALP > ULN versus",
        "ALP <= ULN' (Tamai 2017 Table 2 row label), binarised inline in",
        "model() as alp_high <- (ALP > alp_uln). Tamai 2017 does not print a",
        "numeric ULN, so alp_uln is set to 120 U/L, the same representative",
        "adult cutoff used by the sibling Gupta_2016_lenvatinib.R; this is",
        "the one non-source-derived number in the file and is flagged in the",
        "vignette Assumptions and deviations. Downstream users should either",
        "supply a pre-binarised ALP column (0 / 1 times any positive value",
        "above 120) or edit alp_uln to their laboratory's ULN. Multiplicative",
        "power-form effect on CL/F: 0.852^alp_high, i.e. CL/F 14.8 percent",
        "lower above the ULN (Tamai 2017 Results). Pooled Table 1 ALP is",
        "highly skewed - mean 160.4 U/L (SD 168.8) against a median of only",
        "81.5, range 19.0-1133.0 - and 71 percent of the study 202 HCC",
        "cohort was above the ULN (Tamai 2017 Discussion), which is why this",
        "term matters disproportionately in HCC."
      ),
      source_name = "ALP"
    ),
    CONMED_CYP3A4_IND = list(
      description = "Concomitant CYP3A4 inducer coadministration indicator (1 = coadministered, 0 = not).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant CYP3A4 inducer).",
      notes = paste(
        "Multiplicative power-form effect on CL/F: 1.30^CONMED_CYP3A4_IND",
        "(+30 percent when 1; Tamai 2017 Table 2 theta_INDU = 1.30, RSE",
        "0.534 percent, bootstrap 1.30 (1.23-1.38)). Inducers were drawn",
        "from the strong and moderate CYP3A inducers listed by the FDA",
        "(Tamai 2017 Methods, reference 33); weak inducers were not pooled",
        "in. The exposure window is explicit: coadministration required the",
        "inducer to have been given for at least 10 days on or before the",
        "day of PK assessment, and an inducer stopped within 7 days before",
        "that day still counted as coadministration. Only 16 of 452 subjects",
        "(3.5 percent) were positive (Tamai 2017 Table 1). Tamai 2017",
        "Discussion judges the effect clinically irrelevant because it is",
        "smaller than the 32.6 percent CV interindividual variability in",
        "CL/F."
      ),
      source_name = "INDU"
    ),
    CONMED_CYP3A4_INH = list(
      description = "Concomitant CYP3A inhibitor coadministration indicator (1 = coadministered, 0 = not).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant CYP3A inhibitor).",
      notes = paste(
        "Multiplicative power-form effect on CL/F: 0.922^CONMED_CYP3A4_INH",
        "(-7.8 percent when 1; Tamai 2017 Table 2 theta_INHIB = 0.922, RSE",
        "0.922 percent, bootstrap 0.921 (0.893-0.951)). Inhibitors were",
        "drawn from the strong and moderate CYP3A inhibitors listed by the",
        "FDA (Tamai 2017 Methods, reference 33); weak inhibitors were not",
        "pooled in. Unlike the inducer flag the window is same-day only:",
        "coadministration was defined as the inhibitor being given ON the",
        "day of PK assessment. 34 of 452 subjects (7.5 percent) were",
        "positive (Tamai 2017 Table 1). The paper labels the CYP3A pathway",
        "the predominant route of lenvatinib CYP-mediated metabolism in",
        "human liver microsomes."
      ),
      source_name = "INHIB"
    ),
    DIS_HEALTHY = list(
      description = "Healthy-participant cohort indicator (1 = healthy adult from a phase 1 clinical pharmacology study, 0 = cancer patient).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (cancer patient; the complement is the pooled solid-tumor plus HCC cohort of the 4 phase 1 dose-finding studies and study 202).",
      notes = paste(
        "Multiplicative power-form effect on CL/F: 1.19^DIS_HEALTHY",
        "(+19 percent when 1; Tamai 2017 Table 2 theta_TM = 1.19, RSE 3.26",
        "percent, bootstrap 1.19 (1.11-1.27)). Table 2's abbreviation list",
        "fixes the orientation unambiguously: 'TM, population 0 (cancer",
        "subjects) or 1 (healthy subjects)'. 232 of 452 subjects (51.3",
        "percent) are healthy (Tamai 2017 Table 1 tumor-type row). A",
        "SEPARATE HCC-versus-other-tumor contrast was tested and REJECTED -",
        "Tamai 2017 Results state the HCC effect on CL/F was not",
        "statistically significant once body weight and ALP were in the",
        "model - so this indicator splits healthy from cancer only, and",
        "must not be repurposed as an HCC flag. The +19 percent here is a",
        "little larger than the +15 percent of the sibling",
        "Gupta_2016_lenvatinib.R, which was fitted without study 202."
      ),
      source_name = "TM"
    ),
    FORM_CAPSULE = list(
      description = "Capsule versus tablet formulation indicator (1 = capsule, 0 = tablet).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (tablet; F is anchored at exactly 1 for the tablet arm, which carries no estimated parameter).",
      notes = paste(
        "Relative bioavailability of capsule to tablet is 0.867 (Tamai 2017",
        "Table 2 F1, RSE 1.00 percent, bootstrap 0.867 (0.815-0.908)),",
        "applied to the depot as f(depot) <- (1 - FORM_CAPSULE) +",
        "FORM_CAPSULE * exp(lfdepot). NO interindividual variability is",
        "attached: Tamai 2017 Results state IIV was estimated for every",
        "parameter EXCEPT the two intercompartmental clearances and the",
        "relative bioavailability of capsule to tablet, and the Table 2 IIV",
        "column is '-' on the F1 row. This differs from the sibling",
        "Gupta_2016_lenvatinib.R, which did carry a 30.2 percent CV IIV on",
        "its capsule F. Tamai 2017 does not tabulate the per-arm split of",
        "the 452 subjects between capsule and tablet; the study 202 HCC",
        "cohort and the phase 3 regimen the paper recommends are capsule.",
        "Same orientation as the Gupta 2016 entry in the FORM_CAPSULE",
        "register, i.e. tablet is the F = 1 reference."
      ),
      source_name = "FORM"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "lenvatinib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "lenvatinib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "lenvatinib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "lenvatinib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 452L,
    n_studies = 13L,
    n_observations = paste(
      "8761 lenvatinib plasma concentrations: 5077 from 232 subjects in 8",
      "phase 1 clinical pharmacology studies in healthy adults, 3188 from",
      "155 subjects in 4 phase 1 dose-finding studies in mixed solid tumors",
      "(NCT00121719, NCT00121680, NCT00280397, NCT01268293), and 496 from",
      "the 65 HCC subjects of study 202 (NCT00946153)."
    ),
    age_range = "18.0-85.0 years",
    age_median = "50.0 years (mean 48.7, SD 16.7)",
    weight_range = "42.7-147.0 kg",
    weight_median = "75.1 kg (mean 75.8, SD 17.7)",
    sex_female_pct = 35.8,
    race_ethnicity = c(
      White = 56.0,
      Japanese = 19.9,
      Black = 14.6,
      Other = 7.3,
      Hispanic = 1.3,
      Other_Asian = 0.9
    ),
    disease_state = paste(
      "Pooled cohort: healthy adults (n = 232, 51.3 percent), mixed solid",
      "tumors other than HCC (n = 155, 34.3 percent), and advanced",
      "hepatocellular carcinoma Child-Pugh class A (n = 65, 14.4 percent).",
      "The HCC subjects come from study 202, a multicentre open-label phase",
      "1/2 study whose phase 1 part identified 12 mg once daily in 4-week",
      "continuous cycles as the maximum tolerated dose."
    ),
    dose_range = paste(
      "Study 202 administered 8 to 16 mg/day orally in 4-week cycles; the",
      "phase 2 expansion used 12 mg once daily. Tamai 2017 does not",
      "enumerate the dose levels of the 12 pooled phase 1 studies."
    ),
    regions = "Multiregional; study 202 was conducted in Japan and Korea. Tamai 2017 does not give a per-region breakdown of the pooled analysis set.",
    renal_function = "Creatinine clearance median 103.6 mL/min (mean 103.9, SD 36.3, range 17.0-268.0).",
    hepatic_function = paste(
      "Albumin median 40.0 g/L (range 24.0-52.0); alkaline phosphatase",
      "median 81.5 U/L (range 19.0-1133.0); ALT median 21.0 U/L (range",
      "5.0-660.0); AST median 22.0 U/L (range 8.0-930.0); bilirubin median",
      "10.3 umol/L (range 2.0-101.1). HCC subjects were Child-Pugh class A;",
      "71 percent of them had alkaline phosphatase above the upper limit of",
      "normal."
    ),
    ecog_distribution = "ECOG 0: 91; ECOG 1: 61; ECOG 2: 2; missing: 298 (the 232 healthy subjects account for most of the missing).",
    co_medication = "Concomitant CYP3A4 inducers in 16 of 452 (3.5 percent); concomitant CYP3A inhibitors in 34 of 452 (7.5 percent).",
    notes = paste(
      "Demographics reproduced from Tamai 2017 Table 1 (N = 452). Race",
      "percentages are computed from the Table 1 counts (White 253,",
      "Japanese 90, Black 66, Other 33, Hispanic 6, Other Asian 4), which",
      "sum exactly to 452; sex is 162 female / 290 male. Fitted in NONMEM",
      "7.2 with PDx-Pop 5.0 using FOCE with interaction, and validated by",
      "500 bootstrap replicates whose medians and 2.5-97.5 percentiles are",
      "quoted in the ini() source-trace comments."
    )
  )

  ini({
    # ==================================================================
    # Tamai 2017 Table 2, "Final Model / Population Mean (%RSE)" column.
    # Every value below is the final-model estimate, not the base-model
    # estimate printed in the two columns to its left, and each carries
    # its bootstrap median and 2.5-97.5 percentile interval from the
    # rightmost column as a second, independent confirmation of the
    # transcription.
    #
    # Reference subject for the structural values: 75 kg, cancer patient
    # (TM = 0), no concomitant CYP3A inducer or inhibitor, alkaline
    # phosphatase at or below the ULN, dosed as the tablet formulation.
    # ==================================================================

    # ----- Structural disposition -----
    lcl <- log(6.43); label("Apparent clearance CL/F of the reference subject (L/h)") # Tamai 2017 Table 2 final theta_CL = 6.43 L/h (RSE 2.19%; bootstrap 6.42, 6.07-6.76)
    lvc <- log(47.0); label("Apparent central volume V1/F at 75 kg (L)") # Tamai 2017 Table 2 final theta_V1 = 47.0 L (RSE 4.40%; bootstrap 46.8, 43.9-49.8)
    lvp <- log(31.2); label("Apparent first peripheral volume V2/F at 75 kg (L)") # Tamai 2017 Table 2 final theta_V2 = 31.2 L (RSE 6.76%; bootstrap 31.1, 28.3-33.7)
    lvp2 <- log(34.5); label("Apparent second peripheral volume V3/F at 75 kg (L)") # Tamai 2017 Table 2 final theta_V3 = 34.5 L (RSE 4.14%; bootstrap 34.7, 31.6-37.7)
    lq <- log(3.96); label("Apparent intercompartmental clearance Q1/F between central and peripheral1 at 75 kg (L/h)") # Tamai 2017 Table 2 final theta_Q1 = 3.96 L/h (RSE 3.03%; bootstrap 3.99, 3.57-4.49)
    lq2 <- log(0.726); label("Apparent intercompartmental clearance Q2/F between central and peripheral2 at 75 kg (L/h)") # Tamai 2017 Table 2 final theta_Q2 = 0.726 L/h (RSE 2.91%; bootstrap 0.738, 0.639-0.845)

    # ----- Absorption -----
    # Sequential zero-order release into the depot over D1 followed by
    # first-order absorption at Ka. Tamai 2017 Methods call this
    # "simultaneous first- and 0-order absorption": both processes act on
    # the depot at once while the zero-order window is open. The
    # parameters are named D1 and F1, i.e. the NONMEM duration and
    # bioavailability of dosing COMPARTMENT 1, which is what fixes the
    # dose as entering the depot and only the depot; a parallel split
    # between two input routes would need a fraction parameter that
    # neither Tamai 2017 nor its Table 2 reports. See the vignette
    # Assumptions and deviations.
    lka <- log(1.04); label("First-order absorption rate constant Ka (1/h)") # Tamai 2017 Table 2 final Ka = 1.04 1/h (RSE 6.80%; bootstrap 1.04, 0.933-1.13). The Table 2 row prints the unit as 'L/h', a typographical slip: Ka is a first-order rate constant
    ld1 <- log(1.06); label("Duration of the zero-order release into the depot D1 (h)") # Tamai 2017 Table 2 final D1 = 1.06 h (RSE 5.77%; bootstrap 1.06, 0.987-1.14)
    lfdepot <- log(0.867); label("Relative bioavailability of the capsule versus the tablet reference (unitless)") # Tamai 2017 Table 2 final F1 = 0.867 (RSE 1.00%; bootstrap 0.867, 0.815-0.908)

    # ----- Covariate effects -----
    # The two body-weight powers are ESTIMATED, not fixed at allometric
    # defaults, so they are written bare rather than inside fixed().
    e_wt_cl <- 0.708; label("Power of (WT/75) on CL/F, Q1/F and Q2/F (unitless)") # Tamai 2017 Table 2 final theta_WGT1 = 0.708 (RSE 6.58%; bootstrap 0.711, 0.538-0.886)
    e_wt_vc_vp <- 1.08; label("Power of (WT/75) on V1/F, V2/F and V3/F (unitless)") # Tamai 2017 Table 2 final theta_WGT2 = 1.08 (RSE 5.42%; bootstrap 1.08, 0.876-1.28)

    # Categorical effects are entered in Tamai 2017 as index variables in
    # the power form theta^covariate (Methods, "Categorical covariates
    # were tested and incorporated in the model as index variables"). They
    # are stored here on the log scale so that exp(e * cov) reproduces the
    # source multiplier exactly when cov is 0 or 1.
    e_cyp3a4_ind_cl <- log(1.30); label("Log-effect of concomitant CYP3A4 inducer on CL/F (unitless)") # Tamai 2017 Table 2 final theta_INDU = 1.30 (RSE 0.534%; bootstrap 1.30, 1.23-1.38); Results: CL/F increased by 30%
    e_cyp3a4_inh_cl <- log(0.922); label("Log-effect of concomitant CYP3A inhibitor on CL/F (unitless)") # Tamai 2017 Table 2 final theta_INHIB = 0.922 (RSE 0.922%; bootstrap 0.921, 0.893-0.951); Results: CL/F decreased by 7.8%
    e_dis_healthy_cl <- log(1.19); label("Log-effect of the healthy-participant cohort on CL/F versus cancer patients (unitless)") # Tamai 2017 Table 2 final theta_TM = 1.19 (RSE 3.26%; bootstrap 1.19, 1.11-1.27); Results: healthy subjects had 19% higher CL/F
    e_alp_cl <- log(0.852); label("Log-effect of alkaline phosphatase above the ULN on CL/F (unitless)") # Tamai 2017 Table 2 final theta_ALP = 0.852 (RSE 1.24%; bootstrap 0.855, 0.807-0.910); Results: CL/F decreased by 14.8%

    # ----- Interindividual variability -----
    # Tamai 2017 Table 2 footnote a defines the printed column exactly:
    # 'The %CV for both intersubject and proportional residual variability
    # is an approximation taken as the square root of the variance x 100'.
    # The omega VARIANCE is therefore (%CV/100)^2 directly, and NOT the
    # log-normal conversion log(CV^2 + 1) used by the sibling
    # Gupta_2016_lenvatinib.R. This reading is confirmed independently by
    # the paper's own simulated AUC bounds, which the vignette reproduces
    # to within 0.25% only when exp(omega^2_CL / 2) uses 0.326^2.
    # No IIV is estimated on Q1/F, Q2/F or F1 (Table 2 IIV column is '-'
    # on those three rows, and Results say so in words).
    #
    # CL/F and V1/F are correlated with R = 0.599 (Tamai 2017 Table 2,
    # line below the table: 'Correlation between CL/F and V1/F for final
    # model: R = 0.599'), so the covariance is
    # 0.599 x 0.326 x 0.495 = 0.096661.
    etalcl + etalvc ~ c(
      0.106276,
      0.096661, 0.245025
    ) # Tamai 2017 Table 2 final IIV CL/F 32.6% CV -> 0.326^2; V1/F 49.5% CV -> 0.495^2; correlation R = 0.599
    etalvp ~ 0.389376 # Tamai 2017 Table 2 final IIV V2/F 62.4% CV -> 0.624^2
    etalvp2 ~ 0.176400 # Tamai 2017 Table 2 final IIV V3/F 42.0% CV -> 0.420^2
    etalka ~ 0.216225 # Tamai 2017 Table 2 final IIV Ka 46.5% CV -> 0.465^2
    etald1 ~ 0.467856 # Tamai 2017 Table 2 final IIV D1 68.4% CV -> 0.684^2

    # ----- Residual unexplained variability -----
    # Tamai 2017 uses a THREE-stratum residual model (Methods): a combined
    # additive plus proportional error for time after dose <= 2 h, and
    # separate proportional errors for the phase 1 clinical pharmacology
    # studies and for the cancer-patient studies. Table 2 final estimates
    # are 44.8% CV proportional + 7.35 ng/mL additive for TAD <= 2 h,
    # 17.3% CV proportional for the clinical pharmacology studies, and
    # 30.2% CV proportional for the patient studies. rxode2 cannot switch
    # a residual model on time after dose, so the library model carries
    # the patient-study proportional term - the stratum that applies to
    # the HCC population this paper is about - together with the additive
    # term. The two unused proportional arms are quoted in the vignette
    # Assumptions and deviations. Same simplification, and the same choice
    # of arm, as the sibling Gupta_2016_lenvatinib.R.
    propSd <- 0.302; label("Proportional residual error, cancer-patient studies (fraction)") # Tamai 2017 Table 2 final proportional %CV (patient studies) = 30.2 (RSE 2.00%; bootstrap 30.1, 27.7-32.2)
    addSd <- 7.35; label("Additive residual error, time after dose <= 2 h (ng/mL)") # Tamai 2017 Table 2 final additive (TAD <= 2 h) = 7.35 ng/mL (RSE 16.3%; bootstrap 7.50, 4.62-9.85)
  })

  model({
    # Reference covariate values. 75 kg is written into every Tamai 2017
    # Table 2 covariate equation as '(WGT/75)'.
    ref_wt <- 75

    # Tamai 2017 enters alkaline phosphatase only as an above-/at-or-below-
    # ULN indicator and never prints a numeric ULN. 120 U/L is a
    # representative adult cutoff carried over from the sibling
    # Gupta_2016_lenvatinib.R; it is the single non-source-derived number
    # in this file. See covariateData[[ALP]]$notes and the vignette.
    alp_uln <- 120
    alp_high <- (ALP > alp_uln)

    # Individual parameters. The four categorical / indicator effects and
    # the body-weight power all act multiplicatively on CL/F, reproducing
    # the Table 2 equation
    #   CL/F = theta_CL * (WGT/75)^theta_WGT1 * theta_INDU^INDU
    #          * theta_INHIB^INHIB * theta_TM^TM * theta_ALP^ALP
    cl <- exp(lcl + etalcl) *
      (WT / ref_wt)^e_wt_cl *
      exp(e_cyp3a4_ind_cl * CONMED_CYP3A4_IND) *
      exp(e_cyp3a4_inh_cl * CONMED_CYP3A4_INH) *
      exp(e_dis_healthy_cl * DIS_HEALTHY) *
      exp(e_alp_cl * alp_high)

    # Q1/F and Q2/F share theta_WGT1 with CL/F and carry no IIV.
    q <- exp(lq) * (WT / ref_wt)^e_wt_cl
    q2 <- exp(lq2) * (WT / ref_wt)^e_wt_cl

    # All three volumes share the separately estimated theta_WGT2.
    vc <- exp(lvc + etalvc) * (WT / ref_wt)^e_wt_vc_vp
    vp <- exp(lvp + etalvp) * (WT / ref_wt)^e_wt_vc_vp
    vp2 <- exp(lvp2 + etalvp2) * (WT / ref_wt)^e_wt_vc_vp

    ka <- exp(lka + etalka)
    d1 <- exp(ld1 + etald1)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # Mammillary three-compartment disposition: both peripheral
    # compartments exchange with central. Tamai 2017 Methods list "apparent
    # volume of peripheral compartments (V2/F and V3/F), intercompartmental
    # clearance between V2 and V3 (Q2/F and Q3/F)" in parallel
    # construction, i.e. one intercompartmental clearance per peripheral
    # compartment. The Table 2 abbreviation list renumbers these to Q1 and
    # Q2 and glosses Q2 as "intercompartment clearance between V2 and V3",
    # which would be a catenary chain; that gloss is a carry-over slip from
    # the Methods sentence. The mammillary reading is the one the
    # predecessor model of the same drug and sponsor uses
    # (Gupta_2016_lenvatinib.R) and the one NONMEM's ADVAN12 supplies. See
    # the vignette Assumptions and deviations.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central -
      k12 * central + k21 * peripheral1 -
      k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    # Zero-order release of the dose into the depot over D1. Dose records
    # must carry rate = -2 for rxode2 to use the modelled duration.
    dur(depot) <- d1

    # Bioavailability: exactly 1 on the tablet reference arm, 0.867 on the
    # capsule arm. No IIV on F (Tamai 2017 Table 2 IIV column is '-').
    f(depot) <- (1 - FORM_CAPSULE) + FORM_CAPSULE * exp(lfdepot)

    # Dose in mg and volumes in L give mg/L; x 1000 for the ng/mL scale
    # Tamai 2017 reports concentrations and AUC in.
    Cc <- central / vc * 1000
    Cc ~ prop(propSd) + add(addSd)
  })
}
