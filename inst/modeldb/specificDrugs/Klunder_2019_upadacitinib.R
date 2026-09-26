Klunder_2019_upadacitinib <- function() {
  description <- paste0(
    "Two-compartment population PK model for oral upadacitinib (ABT-494, a ",
    "selective JAK1 inhibitor) pooled across 12 phase I-III trials in healthy ",
    "volunteers and adults with rheumatoid arthritis. The absorption model is ",
    "formulation-dependent: the immediate-release capsule is absorbed ",
    "first-order from a depot with a lag time, while the extended-release ",
    "tablet uses a parallel mixed process in which 74.5% of the absorbed dose ",
    "enters the central compartment by a zero-order input of 3.29 h and the ",
    "remaining 25.5% enters the depot and is absorbed first-order, both arms ",
    "sharing a common lag time and a 76.2% relative bioavailability. ",
    "Statistically significant covariates retained in the final model: ",
    "patient population (rheumatoid arthritis vs healthy) and baseline ",
    "creatinine clearance and bodyweight on CL/F, and bodyweight on Vc/F. ",
    "Intersubject variability on CL/F and Vc/F, and the proportional ",
    "residual-error magnitude, are each estimated separately for the phase I ",
    "studies and for the phase II/III studies. This is the successor to the ",
    "phase I + II immediate-release-only analysis in ",
    "modellib('Klunder_2017_upadacitinib'), and is the parent model whose ",
    "structural parameters were carried forward and fixed in ",
    "modellib('Bhatnagar_2024_upadacitinib')."
  )
  reference <- paste0(
    "Klunder B, Mittapalli RK, Mohamed M-EF, Friedel A, Noertersheuser P, ",
    "Othman AA. Population pharmacokinetics of upadacitinib using the ",
    "immediate-release and extended-release formulations in healthy subjects ",
    "and subjects with rheumatoid arthritis: analyses of phase I-III clinical ",
    "trials. Clin Pharmacokinet. 2019;58(8):1045-1058. ",
    "doi:10.1007/s40262-019-00739-3"
  )
  vignette <- "Klunder_2019_upadacitinib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Baseline total bodyweight.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Power-form effect on both CL/F (exponent 0.132, Klunder 2019 Table 3) ",
        "and Vc/F (exponent 0.804, Klunder 2019 Table 3). Reference bodyweight ",
        "74 kg. Klunder 2019 Methods states only that continuous covariates ",
        "entered 'with a power function centered on the median covariate ",
        "value' and never prints the median itself; the 74 kg divisor is read ",
        "from the NONMEM control stream of the successor analysis that ",
        "inherits this exact model (Bhatnagar 2024 Appendix S1, ",
        "'MU_1 = THETA(1) + THETA(13)*LOG(CRCL/108.70) + THETA(15)*LOG(WTKG/74)' ",
        "and 'MU_2 = THETA(2) + THETA(14)*LOG(WTKG/74)'), corroborated by the ",
        "rendered covariate equation in Bhatnagar 2024 Table S3 footnote b. ",
        "See the vignette Errata. The analysis-dataset bodyweight mean was ",
        "76.4 kg (SD 19.33, range 36-196; Klunder 2019 Table 2)."
      ),
      source_name = "WTKG"
    ),
    CRCL = list(
      description = paste0(
        "Baseline creatinine clearance estimated by the Cockcroft-Gault ",
        "formula, reported as raw mL/min and NOT BSA-normalized."
      ),
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Power-form effect on CL/F (exponent 0.256, Klunder 2019 Table 3). ",
        "Reference CrCL 108.70 mL/min, read from the successor control stream ",
        "'LOG(CRCL/108.70)' (Bhatnagar 2024 Appendix S1) because Klunder 2019 ",
        "does not print its own centering median; see the WT notes and the ",
        "vignette Errata. Like its predecessor ",
        "modellib('Klunder_2017_upadacitinib'), this model uses raw ",
        "Cockcroft-Gault CrCL in mL/min, which is the 'CLCR' / 'CrCL' raw ",
        "branch of the CRCL register entry rather than the BSA-normalized ",
        "mL/min/1.73 m^2 branch -- supplying a BSA-normalized value would bias ",
        "CL/F. The analysis-dataset CrCL mean was 113.7 mL/min (SD 37.92, ",
        "range 30.2-390.9; Klunder 2019 Table 2). Baseline serum bilirubin, ",
        "AST, ALT, age, DAS28-CRP and hsCRP were also tested on CL/F and were ",
        "not retained."
      ),
      source_name = "CRCL"
    ),
    DIS_RA = list(
      description = "Patient-population indicator (1 = subject with rheumatoid arthritis, 0 = healthy volunteer).",
      units = "(binary)",
      type = "binary",
      reference_category = paste0(
        "0 (healthy volunteer). Klunder 2019 Table 3 reports the contrast as ",
        "the 'CL/F ratio of RA patients compared with healthy subjects' = ",
        "0.754, i.e. the healthy volunteer is the reference and the RA patient ",
        "carries the multiplier."
      ),
      notes = paste0(
        "Multiplicative effect on CL/F: RA patients have 0.754 times the ",
        "clearance of healthy subjects (about 25% lower, giving about 33% ",
        "higher AUC; Klunder 2019 Table 3 and Discussion). Encoded on the log ",
        "scale as e_ra_cl = log(0.754) applied as e_ra_cl * DIS_RA. The source ",
        "control stream carries the identical polarity -- 'pref1 = 1; ",
        "IF(RA.EQ.1) pref1 = THETA(12)' with THETA(12) = 0.754 (Bhatnagar 2024 ",
        "Appendix S1) -- so the source column maps onto the canonical ",
        "DIS_RA orientation with no value transformation. The sibling models ",
        "modellib('Klunder_2017_upadacitinib') and ",
        "modellib('Bhatnagar_2024_upadacitinib') carry the same contrast under ",
        "the complementary DIS_HEALTHY column; DIS_RA is used here because ",
        "the patient cohort of this analysis IS the rheumatoid arthritis ",
        "cohort, which is exactly what DIS_RA names, whereas Bhatnagar 2024 ",
        "studies axial spondyloarthritis patients. 3992 of 4170 subjects (96%) ",
        "were RA patients (Klunder 2019 Table 2)."
      ),
      source_name = "RA"
    ),
    FORM_UPA_ER = list(
      description = "Formulation indicator for the dosing record (1 = extended-release tablet, 0 = immediate-release capsule).",
      units = "(binary)",
      type = "binary",
      reference_category = paste0(
        "0 (immediate-release capsule). The immediate-release form is the ",
        "bioavailability reference: Klunder 2019 Table 3 reports the ",
        "extended-release relative bioavailability as 76.2%, so the tabulated ",
        "CL/F, Vc/F, Q/F and Vp/F are all on immediate-release bioavailability."
      ),
      notes = paste0(
        "Selects the entire absorption sub-model, not a single parameter. ",
        "Extended release: Ka = 0.0523 1/h carrying the only absorption random ",
        "effect, a 0.154 h lag on both arms, 74.5% of the absorbed dose ",
        "entering the central compartment over a 3.29 h zero-order input, and ",
        "76.2% relative bioavailability. Immediate release: Ka = 2.77 1/h with ",
        "no random effect, a 0.200 h lag, no zero-order arm and unit relative ",
        "bioavailability (Klunder 2019 Table 3). The source control stream ",
        "implements this as 'IF(FORM.EQ.1) THEN KA = EXP(MM_7); LAG = ",
        "EXP(MM_8); TBIO = 1; INFRAC = 0; ENDIF', i.e. source FORM = 1 is ",
        "immediate release and FORM_UPA_ER = FORM - 1. Doses were 1-48 mg ",
        "immediate release and 7.5-30 mg extended release."
      ),
      source_name = "FORM (1 = immediate release, 2 = extended release; FORM_UPA_ER = FORM - 1)"
    ),
    STUDY_UPA_PHASE1 = list(
      description = "Study-stratum indicator (1 = record from one of the four phase I studies, 0 = record from a phase II, IIb/III or III study).",
      units = "(binary)",
      type = "binary",
      reference_category = paste0(
        "0 (the pooled phase II / IIb-III / III patient studies, which supply ",
        "3982 of the 4170 subjects and are therefore the bulk stratum)."
      ),
      notes = paste0(
        "Switches the interindividual-variability magnitudes on CL/F and Vc/F ",
        "and the proportional residual-error magnitude between the two study ",
        "strata; it does NOT change any typical value. Klunder 2019 Results: ",
        "'A model with different ISV on CL/F and Vc/F between phase I and ",
        "phase II/III studies was found to reduce the OFV by 152 points' and ",
        "'A split of the proportional residual error between phase I and phase ",
        "II/III studies ... improve[d] the fit by a 1108-point reduction in ",
        "the OFV'. The source control stream builds the stratum as a ",
        "'Ph23flag' derived from the study number (Bhatnagar 2024 Appendix ",
        "S1). Phase I contributed 188 subjects; the wider inclusion criteria ",
        "and sparser, less-controlled sampling of the phase II/III trials are ",
        "the paper's stated reason for the larger phase II/III variability ",
        "(Klunder 2019 Discussion)."
      ),
      source_name = "STDY (Ph23flag = 1 for the nine phase II/III study numbers; STUDY_UPA_PHASE1 = 1 - Ph23flag)"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Biological sex indicator (1 = female, 0 = male).",
      units = "(binary)",
      type = "binary",
      notes = paste0(
        "Tested on both CL/F and Vc/F. Sex on Vc/F entered the full model ",
        "during forward inclusion but did not maintain significance ",
        "(p > 0.001) during backward elimination and was removed; sex on CL/F ",
        "was never significant (Klunder 2019 Results). Retained as a ",
        "documented screen only -- the predecessor ",
        "modellib('Klunder_2017_upadacitinib') DID retain sex on both ",
        "parameters, so the disappearance of the effect in the larger ",
        "phase I-III dataset is itself a finding. No point estimate is ",
        "published here, so no value can be carried."
      )
    ),
    RACE_BLACK = list(
      description = "Race indicator (1 = Black, 0 = otherwise).",
      units = "(binary)",
      type = "binary",
      notes = "Screened on CL/F and Vc/F alongside White, Hispanic and Asian indicators; not retained (Klunder 2019 Results and Conclusions). No point estimate published."
    ),
    CONMED_METHOTREXATE = list(
      description = "Concomitant methotrexate use (1 = yes, 0 = no).",
      units = "(binary)",
      type = "binary",
      notes = "Screened on CL/F as a time-varying covariate; not retained. 65% of the analysis population used background methotrexate (Klunder 2019 Table 2). No point estimate published."
    )
  )

  compartmentData <- list(
    depot = list(analyte = "upadacitinib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "upadacitinib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "upadacitinib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 4170,
    n_studies = 12,
    n_observations = 29372,
    age_range = "18-87 years",
    age_mean = "53.9 years (SD 12.73)",
    weight_range = "36-196 kg",
    weight_mean = "76.4 kg (SD 19.33)",
    bmi_mean = "28.5 kg/m^2 (SD 6.60), range 13.3-71.9",
    sex_female_pct = 76,
    race_ethnicity = c(White = 80, Black = 6, Asian = 11, `Multiple races` = 1, Other = 2),
    disease_state = "healthy volunteers (178 subjects, 4%) and adults with moderate to severe rheumatoid arthritis (3992 subjects, 96%)",
    renal_function = "baseline Cockcroft-Gault creatinine clearance mean 113.7 mL/min (SD 37.92), range 30.2-390.9",
    co_medication = "background methotrexate in 65%; pH-modifying agents in 40%; moderate CYP3A inhibitors in 3%",
    dose_range = "1-48 mg immediate-release capsule (single and twice-daily) and 7.5-30 mg extended-release tablet (once daily)",
    regions = "global, including a regional phase IIb/III study in Japanese subjects",
    notes = paste0(
      "Four phase I studies (healthy volunteers plus 14 subjects with mild to ",
      "moderate RA), two phase II studies, one regional Japanese phase IIb/III ",
      "study and five global phase III studies. Baseline demographics are ",
      "Klunder 2019 Table 2, 'All subjects' column; the per-phase columns of ",
      "that table give the phase I subgroup separately (n = 188, mean age 36.0 ",
      "years, mean weight 75.4 kg, 87% male)."
    )
  )

  ini({
    # ======================================================================
    # All values are Klunder 2019 Table 3, column 'Population analysis /
    # Estimate (%RSE)'. Every parameter in that table was ESTIMATED -- the
    # table reports a %RSE and a bootstrap 95% CI for each -- so nothing here
    # is wrapped in fixed(). (The successor axSpA analysis in
    # modellib('Bhatnagar_2024_upadacitinib') fixes most of these same values;
    # that is a property of the successor, not of this parent model.)
    #
    # SOURCE-TRACE CONFIRMATION. Table 3 is reported on IMMEDIATE-RELEASE
    # bioavailability. The abstract instead quotes the extended-release
    # numbers for healthy volunteers, and the two agree exactly through the
    # 76.2% relative bioavailability:
    #   CL/F(ER)  = 40.9 / 0.762        = 53.7 L/h   (abstract: 53.7 L/h)
    #   Vss/F(ER) = (156 + 68.0) / 0.762 = 294 L     (abstract: 294 L)
    # ======================================================================

    # ---------------- Disposition ------------------------------------------
    lcl <- log(40.9)
    label("Apparent oral clearance CL/F in a healthy volunteer at the reference covariates (L/h)")
    # Table 3 'CL/F (L/h)' = 40.9 (%RSE 1.6); bootstrap median 41.3, 95% CI 39.6-42.5.

    lvc <- log(156)
    label("Apparent central volume of distribution Vc/F at the reference bodyweight (L)")
    # Table 3 'Vc/F (L)' = 156 (%RSE 1.7); bootstrap median 156, 95% CI 150-161.

    lq <- log(3.22)
    label("Apparent intercompartmental clearance Q/F (L/h)")
    # Table 3 'Q/F (L/h)' = 3.22 (%RSE 5.8); bootstrap median 3.22, 95% CI 2.86-3.63.

    lvp <- log(68.0)
    label("Apparent peripheral volume of distribution Vp/F (L)")
    # Table 3 'Vp/F (L)' = 68.0 (%RSE 7.2); bootstrap median 67.4, 95% CI 59.7-78.3.

    # ---------------- Extended-release absorption ---------------------------
    lka_er <- log(0.0523)
    label("First-order absorption rate constant Ka of the extended-release tablet (1/h)")
    # Table 3 'Extended-release Ka (1/h)' = 0.0523 (%RSE 6.0); bootstrap median 0.0523, 95% CI 0.0460-0.0590.

    ltlag_er <- log(0.154)
    label("Absorption lag time of the extended-release tablet (h)")
    # Table 3 'Extended-release absorption lag time (h)' = 0.154 (%RSE 7.7). Applied to BOTH extended-release arms: the source control stream sets ALAG1 = ALAG2 = LAG (Bhatnagar 2024 Appendix S1).

    logitffo <- log(0.255 / 0.745)
    label("Logit of the fraction of the absorbed extended-release dose taking the FIRST-ORDER depot route (unitless)")
    # Complement of Table 3 'Fraction of extended-release dose absorbed through zero-order process (%)' = 74.5 (%RSE 1.7), so the first-order share is 1 - 0.745 = 0.255 and logit(0.255) = -1.0721. The source estimates this on the logit scale: the control stream sets INFRAC = EXP(THETA(5))/(1+EXP(THETA(5))) with THETA(5) = 1.07, and expit(1.07) = 0.7446 (Bhatnagar 2024 Appendix S1). logit(1-x) = -logit(x), so carrying the first-order share is the same parameter with the sign flipped.

    ld2 <- log(3.29)
    label("Duration of the zero-order absorption input into the central compartment (h)")
    # Table 3 'Zero-order infusion duration (h)' = 3.29 (%RSE 1.7). Enters the source control stream as D2, i.e. the duration on the CENTRAL compartment; the supplement names the row 'Zero-Order Infusion duration D2 (h)'.

    lfrel_er <- log(0.762)
    label("Bioavailability of the extended-release tablet relative to the immediate-release capsule (fraction)")
    # Table 3 'Bioavailability of the extended-release formulation relative to the immediate-release formulation (%)' = 76.2 (%RSE 1.4).

    # ---------------- Immediate-release absorption --------------------------
    lka_ir <- log(2.77)
    label("First-order absorption rate constant Ka of the immediate-release capsule (1/h)")
    # Table 3 'Immediate-release Ka (1/h)' = 2.77 (%RSE 7.4).

    ltlag_ir <- log(0.200)
    label("Absorption lag time of the immediate-release capsule (h)")
    # Table 3 'Immediate-release absorption lag time (h)' = 0.200 (%RSE 3.9).

    # ---------------- Covariate effects -------------------------------------
    # Control-stream forms, on the log scale before exponentiation:
    #   MU_1 = THETA(1) + THETA(13)*LOG(CRCL/108.70) + THETA(15)*LOG(WTKG/74)
    #   MU_2 = THETA(2) + THETA(14)*LOG(WTKG/74)
    #   pref1 = 1; IF(RA.EQ.1) pref1 = THETA(12); CL = EXP(MU_1 + ETA(1))*pref1
    e_crcl_cl <- 0.256
    label("Power exponent: baseline creatinine clearance on CL/F (unitless)")
    # Table 3 'Covariate exponent of creatinine clearance on CL/F' = 0.256 (%RSE 10.0).

    e_wt_cl <- 0.132
    label("Power exponent: bodyweight on CL/F (unitless)")
    # Table 3 'Covariate exponent of weight on CL/F' = 0.132 (%RSE 28.7) -- the least precisely estimated fixed effect in the model.

    e_wt_vc <- 0.804
    label("Power exponent: bodyweight on Vc/F (unitless)")
    # Table 3 'Covariate exponent of weight on Vc/F' = 0.804 (%RSE 8.0).

    e_ra_cl <- log(0.754)
    label("Log ratio of CL/F in rheumatoid arthritis patients relative to healthy volunteers (unitless)")
    # Table 3 'CL/F ratio of RA patients compared with healthy subjects' = 0.754 (%RSE 1.7). Discussion: 'RA subjects were estimated to have 25% lower upadacitinib clearance (leading to 33% higher estimated upadacitinib AUC) compared with healthy subjects.'

    # ---------------- Interindividual variability ---------------------------
    # Variances on the log scale. Table 3 reports each ISV as a percentage
    # under the footnote '%ISV was calculated as SQRT(omega^2) x 100', so the
    # variance is the square of the tabulated percentage divided by 100.
    # CL/F and Vc/F each carry a phase I and a phase II/III magnitude; the
    # STUDY_UPA_PHASE1 indicator selects which applies to a given subject.
    # No CL/F-Vc/F covariance is reported for this model.
    etalcl_ph1 ~ 0.042025
    # Table 3 'ISV on CL/F in phase I (%)' = 20.5; 0.205^2 = 0.042025

    etalcl_ph23 ~ 0.133225
    # Table 3 'ISV on CL/F in phase II/III (%)' = 36.5; 0.365^2 = 0.133225

    etalvc_ph1 ~ 0.059536
    # Table 3 'ISV on Vc/F in phase I (%)' = 24.4; 0.244^2 = 0.059536

    etalvc_ph23 ~ 0.280900
    # Table 3 'ISV on Vc/F in phase II/III (%)' = 53.0; 0.530^2 = 0.2809

    etalka_er ~ 0.446224
    # Table 3 'ISV on extended-release Ka (%)' = 66.8; 0.668^2 = 0.446224. Applies to the EXTENDED-RELEASE arm only -- Klunder 2019 Results states that 'the estimation of ISV variability on the Ka of the IR formulation was not numerically feasible with the current dataset'.

    # ---------------- Residual error ----------------------------------------
    # Source $ERROR block: Y = IPRED + IPRED*EPS(1) + EPS(2) with a diagonal
    # $SIGMA, i.e. a combined proportional-plus-additive model in which the
    # two variances add -- nlmixr2's default add() + prop() parameterization.
    # Klunder 2019 Eq. 2 states the same structure.
    propSdPhase1 <- 0.344
    label("Proportional residual error SD, phase I stratum (fraction)")
    # Table 3 'Proportional error SD in phase I' = 0.344 (%RSE 23.9).

    propSdPhase23 <- 0.543
    label("Proportional residual error SD, phase II/III stratum (fraction)")
    # Table 3 'Proportional error SD in phase II/III' = 0.543 (%RSE 14.0).

    addSd <- 0.0858
    label("Additive residual error SD (ng/mL)")
    # Table 3 'Additive error SD (ng/mL)' = 0.0858 (%RSE 54.5). Shared across both study strata.
  })

  model({
    # ----------------------------------------------------------------------
    # 1. Reference covariate values.
    #
    #    Klunder 2019 Methods says only that continuous covariates entered
    #    "with a power function centered on the median covariate value" and
    #    never prints either median. Both divisors below are read from the
    #    NONMEM control stream of the successor analysis that inherits this
    #    model unchanged (Bhatnagar 2024 Appendix S1), and are corroborated
    #    by the rendered covariate equation in Bhatnagar 2024 Table S3
    #    footnote b. See the vignette Errata.
    # ----------------------------------------------------------------------
    ref_wt <- 74 # kg
    ref_crcl <- 108.70 # mL/min

    # ----------------------------------------------------------------------
    # 2. Formulation-gated absorption.
    #
    #    Reproduces the control stream's
    #      IF(FORM.EQ.1) THEN KA = EXP(MM_7); LAG = EXP(MM_8); TBIO = 1;
    #                         INFRAC = 0; ENDIF
    #    with indicator arithmetic rather than a branch, so the
    #    extended-release terms -- including the Ka random effect, which the
    #    immediate-release branch does not carry -- switch off wholesale when
    #    FORM_UPA_ER = 0.
    #
    #    ffo is the share of the absorbed dose taking the FIRST-ORDER depot
    #    route. For the immediate release capsule there is no zero-order arm,
    #    so ffo collapses to 1 and the whole dose goes through the depot;
    #    exp(0) = 1 likewise recovers unit relative bioavailability.
    # ----------------------------------------------------------------------
    ka <- exp((lka_er + etalka_er) * FORM_UPA_ER + lka_ir * (1 - FORM_UPA_ER))
    tlag <- exp(ltlag_er * FORM_UPA_ER + ltlag_ir * (1 - FORM_UPA_ER))
    frel <- exp(lfrel_er * FORM_UPA_ER)
    ffo <- expit(logitffo) * FORM_UPA_ER + (1 - FORM_UPA_ER)

    # ----------------------------------------------------------------------
    # 3. Phase-stratified interindividual variability.
    #
    #    Only the VARIANCE differs between the two study strata -- the
    #    typical values of CL/F and Vc/F are common to both -- so a single
    #    structural mean is gated against a pair of random effects.
    # ----------------------------------------------------------------------
    etalcl_eff <- etalcl_ph1 * STUDY_UPA_PHASE1 + etalcl_ph23 * (1 - STUDY_UPA_PHASE1)
    etalvc_eff <- etalvc_ph1 * STUDY_UPA_PHASE1 + etalvc_ph23 * (1 - STUDY_UPA_PHASE1)

    # ----------------------------------------------------------------------
    # 4. Individual parameters.
    # ----------------------------------------------------------------------
    cl <- exp(
      lcl + etalcl_eff +
        e_crcl_cl * log(CRCL / ref_crcl) +
        e_wt_cl * log(WT / ref_wt) +
        e_ra_cl * DIS_RA
    )
    vc <- exp(lvc + etalvc_eff + e_wt_vc * log(WT / ref_wt))
    q <- exp(lq)
    vp <- exp(lvp)
    d2 <- exp(ld2)

    # ----------------------------------------------------------------------
    # 5. Micro-constants for the explicit two-compartment form (the source
    #    uses NONMEM ADVAN4: K = CL/V2, K23 = Q/V2, K32 = Q/V3).
    # ----------------------------------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ----------------------------------------------------------------------
    # 6. ODE system: two-compartment disposition with a first-order depot arm
    #    and a parallel zero-order input straight into the central
    #    compartment.
    #
    #    DOSING. Each EXTENDED-RELEASE administration is TWO dose records in
    #    the event table at the same time and the same nominal amount:
    #      - one to cmt = "depot"                  (first-order arm)
    #      - one to cmt = "central" with rate = -2 (zero-order arm; rate = -2
    #        is what makes rxode2 honour the modelled dur(central) instead of
    #        delivering a bolus)
    #    An IMMEDIATE-RELEASE administration needs only the depot record,
    #    because ffo = 1 makes f(central) zero.
    # ----------------------------------------------------------------------
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - (kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # ----------------------------------------------------------------------
    # 7. Dose partitioning, zero-order duration and lag time.
    #    Control stream: F1 = TBIO*(1-INFRAC), F2 = TBIO*INFRAC, D2 = INDUR,
    #    ALAG1 = ALAG2 = LAG. Here ffo = 1 - INFRAC.
    # ----------------------------------------------------------------------
    f(depot) <- frel * ffo
    f(central) <- frel * (1 - ffo)
    dur(central) <- d2
    alag(depot) <- tlag
    alag(central) <- tlag

    # ----------------------------------------------------------------------
    # 8. Observation. Doses are in mg and vc is in L, so central/vc is mg/L;
    #    the source reports plasma concentrations in ng/mL, hence the factor
    #    of 1000. The additive residual SD of 0.0858 is likewise ng/mL.
    #
    #    The proportional residual magnitude is selected per record by the
    #    same study-stratum indicator that selects the ISV magnitudes.
    # ----------------------------------------------------------------------
    propSd <- propSdPhase1 * STUDY_UPA_PHASE1 + propSdPhase23 * (1 - STUDY_UPA_PHASE1)

    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
