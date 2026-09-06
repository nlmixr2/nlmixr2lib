Bhatnagar_2024_upadacitinib <- function() {
  description <- paste0(
    "Two-compartment population PK model for oral upadacitinib (ABT-494, a ",
    "selective JAK1 inhibitor) in adults with axial spondyloarthritis ",
    "(axSpA: ankylosing spondylitis and non-radiographic axSpA) from the ",
    "SELECT-AXIS 1 and SELECT-AXIS 2 trials. Absorption of the ",
    "extended-release tablet is a parallel mixed process: a fraction of the ",
    "dose enters the central compartment by a zero-order input of fixed ",
    "duration and the remainder enters a depot compartment absorbed ",
    "first-order, both delayed by a common lag time. The model carries a ",
    "formulation switch so the same file also reproduces the ",
    "immediate-release capsule arm of the parent analysis (first-order ",
    "absorption only, a different Ka and lag time, and unit relative ",
    "bioavailability). Every structural parameter except the apparent ",
    "central volume was fixed from the previously published upadacitinib ",
    "population PK model built on 4170 healthy volunteers and patients with ",
    "rheumatoid arthritis; the central volume, the interindividual ",
    "variability and the residual error were re-estimated on the ",
    "SELECT-AXIS 1 ankylosing spondylitis data and then held fixed for ",
    "SELECT-AXIS 2. Retained covariates are creatinine clearance and body ",
    "weight on apparent oral clearance, body weight on apparent central ",
    "volume, and a patient-versus-healthy-volunteer clearance ratio. ",
    "Companion exposure-response models from the same paper cover the ASAS20 ",
    "and ASAS40 week-14 response endpoints in each patient population."
  )
  reference <- paste(
    "Bhatnagar S, Eckert D, Stodtmann S, Song I-H, Wung P, Liu W,",
    "Mohamed M-EF. Population pharmacokinetics and exposure-response",
    "analyses for efficacy and safety of upadacitinib in patients with",
    "axial spondyloarthritis. Clin Transl Sci. 2024;17(2):e13733.",
    "doi:10.1111/cts.13733.",
    "Parameter values are Table S3 and the SELECT-AXIS-2 NONMEM control",
    "stream in Appendix S1 of the Supporting Information. The fixed",
    "structural parameters and covariate coefficients originate in the",
    "upstream pooled healthy-volunteer / rheumatoid-arthritis population PK",
    "model of Klunder B, Mittapalli RK, Mohamed M-EF, et al. (Bhatnagar 2024",
    "reference 17), which is not itself in nlmixr2lib; the related earlier",
    "phase I/II analysis is available as modellib('Klunder_2017_upadacitinib').",
    sep = " "
  )
  vignette <- "Bhatnagar_2024_upadacitinib"
  units    <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot       = list(analyte = "upadacitinib", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "upadacitinib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "upadacitinib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight at baseline.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effects on BOTH apparent oral clearance (exponent 0.132) and apparent central volume (exponent 0.804), each normalized to a 74 kg reference. The reference weight and both exponents are fixed from the upstream pooled healthy-volunteer / rheumatoid-arthritis model (Bhatnagar 2024 reference 17) and appear as THETA(15) and THETA(14) in the Appendix S1 control streams. Observed range in the axSpA subjects with PK sampling was 41.5-144 kg (Table S1), so the exponents are exercised well outside the reference. Note that the Table S3 footnote b equation prints 0.123 for the clearance exponent and 0.864 for the volume exponent, which are neither the tabulated nor the control-stream values -- see the vignette Errata.",
      source_name        = "WTKG"
    ),
    CRCL = list(
      description        = "Baseline creatinine clearance (raw mL/min, NOT BSA-normalized).",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effect on apparent oral clearance (exponent 0.256) normalized to a 108.70 mL/min reference; both are fixed from the upstream model and appear as THETA(13) and the literal 108.70 divisor in the Appendix S1 control streams. Bhatnagar 2024 does not restate the estimating equation, but the reference value and the units carry over from the upstream upadacitinib analyses, which used raw Cockcroft-Gault creatinine clearance in mL/min rather than a BSA-normalized eGFR; the companion registered model Klunder_2017_upadacitinib.R uses the same raw form with a 107 mL/min reference. Do NOT supply a mL/min/1.73 m^2 value here. Bhatnagar 2024 reports that upadacitinib exposures are comparable between patients with mild or moderate renal impairment and those with normal renal function despite the statistically significant covariate.",
      source_name        = "CRCL"
    ),
    DIS_HEALTHY = list(
      description        = "Healthy-participant cohort indicator (1 = healthy volunteer, 0 = patient).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (healthy volunteer). The structural lcl is the healthy-volunteer clearance; patients carry the additional multiplicative ratio 0.754 (Table S3, 'CL/F Ratio of Diseased Patients Compared to Healthy Patients'), i.e. clearance is about 25% lower in patients, as stated in the Bhatnagar 2024 Results.",
      notes              = "Time-fixed per subject. EVERY subject in the SELECT-AXIS analyses is a patient, so the axSpA fit is entirely at DIS_HEALTHY = 0: both Appendix S1 control streams hard-code RA = 1 and therefore apply THETA(12) = 0.754 to every record. The indicator is retained here because the fixed structural clearance is expressed on the healthy-volunteer scale inherited from the upstream pooled model, so reproducing any published axSpA exposure requires DIS_HEALTHY = 0. Same canonical orientation as Klunder_2017_upadacitinib.R.",
      source_name        = "RA (1 = patient in the control stream; DIS_HEALTHY = 1 - RA)"
    ),
    FORM_UPA_ER = list(
      description        = "Upadacitinib formulation indicator (1 = extended-release tablet, 0 = immediate-release capsule).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (immediate-release capsule). Table S3 footnote a states that the tabulated CL/F, Vc/F, Q/F and Vp/F 'are for the immediate-release formulation (based on immediate-release bioavailability)', so the immediate-release arm is the bioavailability reference and the extended-release arm carries the 0.762 relative-bioavailability factor.",
      notes              = "Per dose occasion. Selects the whole absorption branch of the model: the extended-release arm uses Ka = 0.0523 1/h with its own interindividual variability, a 0.154 h lag time, 74.5% of the absorbed dose entering the central compartment by a 3.29 h zero-order input, and 76.2% relative bioavailability; the immediate-release arm uses Ka = 2.77 1/h with no interindividual variability, a 0.200 h lag time, no zero-order arm and unit relative bioavailability. Every axSpA subject in SELECT-AXIS 1 and 2 received the 15 mg once-daily EXTENDED-RELEASE tablet, so FORM_UPA_ER = 1 reproduces this paper; the SELECT-AXIS-2 control stream hard-codes FORM = 2 (extended release). The immediate-release branch is retained because it is fully parameterized in the Appendix S1 control streams and in Table S3, and because it defines the bioavailability reference the clearance and volume estimates are expressed on. With FORM_UPA_ER = 0 the zero-order fraction collapses to zero, so an immediate-release dose needs only the depot dose record.",
      source_name        = "FORM (1 = immediate release, 2 = extended release; FORM_UPA_ER = FORM - 1)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 244L,
    n_studies      = 3L,
    age_range      = "19-82 years (subjects with PK sampling, Table S1)",
    age_median     = "44 years (subjects with PK sampling, Table S1; mean 45, SD 13)",
    weight_range   = "41.5-144 kg (subjects with PK sampling, Table S1)",
    weight_median  = "79.0 kg (subjects with PK sampling, Table S1; mean 81.2, SD 20.0). The covariate reference weight of 74 kg is inherited from the upstream pooled model, not from this cohort.",
    sex_female_pct = 45,
    race_ethnicity = "Not tabulated. Bhatnagar 2024 Figure 4 compares model-estimated average concentrations between Asian and non-Asian patients and concludes race has no clinically meaningful effect, but reports no group sizes or percentages.",
    disease_state  = "Axial spondyloarthritis. The population PK datasets comprised 173 patients with ankylosing spondylitis (SELECT-AXIS 1, bDMARD-naive; and SELECT-AXIS 2 study 1, bDMARD inadequate responders) and 71 patients with non-radiographic axial spondyloarthritis (SELECT-AXIS 2 study 2).",
    dose_range     = "Upadacitinib 15 mg once daily, extended-release tablet, versus placebo (1:1 randomization) through the week-14 primary end point in all three studies.",
    regions        = "Multinational phase II/III programme; not broken out in Bhatnagar 2024.",
    renal_function = "Creatinine clearance was a retained covariate on apparent oral clearance; Bhatnagar 2024 reports that plasma exposures are comparable between patients with mild or moderate renal impairment and patients with normal renal function. The cohort renal-function range is not tabulated.",
    notes          = "Clinical trial registrations NCT03178487 (SELECT-AXIS 1) and NCT04169373 (SELECT-AXIS 2 studies 1 and 2). PK samples were collected in all SELECT-AXIS 1 patients and in about 30% of SELECT-AXIS 2 patients at weeks 2, 8, 12 and 14; assay lower limit of quantitation 0.05 ng/mL, with below-limit samples imputed at LLOQ/2 (M5). Table S1 summarizes all 730 randomized subjects, of whom 295 had PK sampling scheduled; the population PK analysis datasets contained the 244 patients with evaluable concentrations quoted here. Fitted in NONMEM with ADVAN4 and FOCE-I. The only difference between the subjects with and without PK sampling was that mean C-reactive protein was about 40% lower in the sampled group; Bhatnagar 2024 notes that CRP has not affected upadacitinib PK in other rheumatological populations."
  )

  ini({
    # ======================================================================
    # All values are Bhatnagar 2024 Table S3 with the corresponding THETA of
    # the SELECT-AXIS-2 control stream in Appendix S1 given for each. The
    # SELECT-AXIS-2 run is the final model: it is a MAXEVAL=0 evaluation, so
    # every value below is what the published axSpA model actually used.
    #
    # fixed() marks exactly the rows Table S3 flags "(FIX)". The apparent
    # central volume, the three interindividual variances and the two
    # residual-error terms are the estimated rows (Table S3 reports a 95%
    # confidence interval for Vc/F and percent relative standard errors for
    # the rest); they were estimated on the SELECT-AXIS 1 ankylosing
    # spondylitis data and, per Table S3 footnote c, carried forward without
    # re-estimation for SELECT-AXIS 2.
    #
    # SOURCE-TRACE CONFIRMATION. For a typical patient (DIS_HEALTHY = 0) the
    # model gives CL/F = 40.9 * 0.754 = 30.8 L/h, so on 15 mg once daily of
    # the extended-release tablet
    #   Cavg = 0.762 * 15 mg * 1000 / (30.8 L/h * 24 h) = 15.4 ng/mL,
    # against the Table 1 medians of 14.5 ng/mL (AS) and 14.8 ng/mL
    # (nr-axSpA). The residual gap is the cohort's above-reference body
    # weight and creatinine clearance, both of which raise clearance.
    # ======================================================================

    # ---------------- Disposition ------------------------------------------
    lcl <- fixed(log(40.9))
    label("Apparent oral clearance CL/F in a healthy volunteer at the reference covariates (L/h)")
    # Table S3 "CL/F (L/h) 40.9 (FIX)"; control-stream THETA(1) = 3.71, exp(3.71) = 40.85. Footnote a: expressed on immediate-release bioavailability.

    lvc <- log(171)
    label("Apparent central volume of distribution Vc/F at the reference body weight (L)")
    # Table S3 "Vc/F (L) 171", 95% CI 128-227; control-stream THETA(2) = 5.14, exp(5.14) = 170.7. The ONLY re-estimated structural parameter: Bhatnagar 2024 Results states re-estimating it "improved the stability of the model", and 171 L is "similar to estimate of 156 L obtained from the previously developed RA model".

    lq <- fixed(log(3.22))
    label("Apparent intercompartmental clearance Q/F (L/h)")
    # Table S3 "Q/F (L/h) 3.22 (FIX)"; control-stream THETA(10) = 1.17, exp(1.17) = 3.222.

    lvp <- fixed(log(68.0))
    label("Apparent peripheral volume of distribution Vp/F (L)")
    # Table S3 "Vp/F (L) 68.0 (FIX)"; control-stream THETA(11) = 4.22, exp(4.22) = 68.03.

    # ---------------- Extended-release absorption ---------------------------
    lka_er <- fixed(log(0.0523))
    label("First-order absorption rate constant Ka of the extended-release tablet (1/h)")
    # Table S3 "Extended-Release KA (1/h) 0.0523 (FIX)"; control-stream THETA(3) = -2.95, exp(-2.95) = 0.05234.

    ltlag_er <- fixed(log(0.154))
    label("Absorption lag time of the extended-release tablet (h)")
    # Table S3 "Extended-Release Lag time (h) 0.154 (FIX)"; control-stream THETA(4) = -1.87, exp(-1.87) = 0.1541. Applied to BOTH absorption arms (ALAG1 = ALAG2 = LAG).

    logitalpha <- fixed(qlogis(0.745))
    label("Logit of the fraction of the absorbed extended-release dose entering by the zero-order process (unitless)")
    # Table S3 "Fraction of Extended-Release Dose Absorbed through Zero-Order Process (%) 74.5 (FIX)"; control-stream THETA(5) = 1.07 entering as INFRAC = exp(THETA5)/(1+exp(THETA5)) = 0.7446, i.e. the source itself parameterizes this fraction on the logit scale.

    ld2 <- fixed(log(3.29))
    label("Duration of the zero-order absorption input into the central compartment (h)")
    # Table S3 "Zero-Order Infusion Duration (h) 3.29 (FIX)"; control-stream THETA(6) = 1.19, exp(1.19) = 3.287, entering as D2.

    lfrel_er <- fixed(log(0.762))
    label("Bioavailability of the extended-release tablet relative to the immediate-release capsule (fraction)")
    # Table S3 "Bioavailability of the Extended-Release Formulation Relative to the Immediate-Release Formulation (%) 76.2 (FIX)"; control-stream THETA(9) = -0.272, exp(-0.272) = 0.7618, entering as TBIO. This is F_rel of Bhatnagar 2024 Equation 1.

    # ---------------- Immediate-release absorption --------------------------
    lka_ir <- fixed(log(2.77))
    label("First-order absorption rate constant Ka of the immediate-release capsule (1/h)")
    # Table S3 "Immediate-Release KA (1/h) 2.77 (FIX)"; control-stream THETA(7) = 1.02, exp(1.02) = 2.773.

    ltlag_ir <- fixed(log(0.200))
    label("Absorption lag time of the immediate-release capsule (h)")
    # Control-stream THETA(8) = -1.61, exp(-1.61) = 0.1999. Table S3 prints "Immediate-Release Lag time (h) 2.00 (FIX)", which is a factor of ten away from its own control stream; the control-stream value is used. See the vignette Errata.

    # ---------------- Covariate effects -------------------------------------
    # Control-stream forms, both on the log scale before exponentiation:
    #   MU_1 = THETA(1) + THETA(13)*LOG(CRCL/108.70) + THETA(15)*LOG(WTKG/74)
    #   MU_2 = THETA(2) + THETA(14)*LOG(WTKG/74)
    #   CL   = EXP(MU_1 + ETA(1)) * THETA(12)^(patient)
    e_crcl_cl <- fixed(0.256)
    label("Power exponent: creatinine clearance on CL/F (unitless)")
    # Table S3 "Covariate Exponent of Creatinine Clearance on CL/F 0.256 (FIX)"; control-stream THETA(13) = 2.56E-01.

    e_wt_cl <- fixed(0.132)
    label("Power exponent: body weight on CL/F (unitless)")
    # Table S3 "Covariate Exponent of Weight on CL/F 0.132 (FIX)"; control-stream THETA(15) = 1.32E-01.

    e_wt_vc <- fixed(0.804)
    label("Power exponent: body weight on Vc/F (unitless)")
    # Table S3 "Covariate Exponent of Weight on Vc/F 0.804 (FIX)"; control-stream THETA(14) = 8.04E-01.

    e_patient_cl <- fixed(log(0.754))
    label("Log ratio of CL/F in patients relative to healthy volunteers (unitless)")
    # Table S3 "CL/F Ratio of Diseased Patients Compared to Healthy Patients 0.754 (FIX)"; control-stream THETA(12) = 7.54E-01 applied as pref1. Bhatnagar 2024 Results: "The clearance of upadacitinib was estimated to be ~25% lower in patients compared to healthy volunteers."

    # ---------------- Interindividual variability ---------------------------
    # Variances on the log scale, from the SELECT-AXIS-2 control-stream
    # $OMEGA blocks. Table S3 reports only the diagonal, as
    # "%IIV was calculated as SQRT(omega^2) x 100":
    #   sqrt(0.111) = 33%, sqrt(0.593) = 77%, sqrt(0.643) = 80%,
    # matching the tabulated 33, 77 and 80. The CL/F-Vc/F COVARIANCE of
    # -0.160 (correlation -0.62) is reported ONLY in the control stream; it
    # is not in Table S3.
    etalcl + etalvc ~ c(0.111, -0.160, 0.593)
    # Control-stream $OMEGA BLOCK(2): 1.11E-01 / -1.60E-01 5.93E-01. Table S3 "IIV on CL/F (%) 33 (56)" and "IIV on Vc/F (%) 77 (70)", parenthetical values being percent relative standard error.

    etalka_er ~ 0.643
    # Control-stream second $OMEGA: 6.43E-01. Table S3 "IIV on Extended-Release KA (%) 80 (57)". Applies to the EXTENDED-RELEASE absorption arm only -- the control stream writes the immediate-release branch as EXP(MM_7) with no eta.

    # ---------------- Residual error ----------------------------------------
    # Control-stream $ERROR: Y = IPRED + IPRED*EPS(1) + EPS(2), i.e. a
    # combined proportional-plus-additive model on the ng/mL scale.
    propSd <- 0.559
    label("Proportional residual error (fraction)")
    # Table S3 "Proportional Error in Phase 3 SD 0.559 (26)"; control-stream $SIGMA 3.12E-01, sqrt(0.312) = 0.5586.

    addSd <- 0.00244
    label("Additive residual error (ng/mL)")
    # Table S3 "Additive Error SD (ng/mL) 0.00244 (110)"; control-stream $SIGMA 5.95E-06, sqrt(5.95e-06) = 0.002439.
  })

  model({
    # ----------------------------------------------------------------------
    # 1. Reference covariate values, both taken verbatim from the Appendix S1
    #    control streams (the LOG(CRCL/108.70) and LOG(WTKG/74) divisors).
    # ----------------------------------------------------------------------
    ref_crcl <- 108.70  # mL/min
    ref_wt   <- 74      # kg

    # ----------------------------------------------------------------------
    # 2. Individual parameters.
    #
    #    The formulation switch reproduces the control stream's
    #      IF(FORM.EQ.1) THEN KA = EXP(MM_7); LAG = EXP(MM_8); TBIO = 1;
    #                         INFRAC = 0; ENDIF
    #    block using indicator arithmetic rather than a branch, so the
    #    extended-release terms (including the Ka random effect, which the
    #    immediate-release branch does not carry) are switched off wholesale
    #    when FORM_UPA_ER = 0.
    # ----------------------------------------------------------------------
    ka   <- exp((lka_er + etalka_er) * FORM_UPA_ER + lka_ir * (1 - FORM_UPA_ER))
    tlag <- exp(ltlag_er * FORM_UPA_ER + ltlag_ir * (1 - FORM_UPA_ER))

    # exp(0) = 1 recovers the immediate-release reference bioavailability,
    # and multiplying the zero-order fraction by the indicator collapses the
    # immediate-release dose onto the first-order arm alone.
    ftot  <- exp(lfrel_er * FORM_UPA_ER)
    alpha <- expit(logitalpha) * FORM_UPA_ER

    cl <- exp(lcl + etalcl +
                e_crcl_cl * log(CRCL / ref_crcl) +
                e_wt_cl * log(WT / ref_wt) +
                e_patient_cl * (1 - DIS_HEALTHY))
    vc <- exp(lvc + etalvc + e_wt_vc * log(WT / ref_wt))
    vp <- exp(lvp)
    q  <- exp(lq)
    d2 <- exp(ld2)

    # ----------------------------------------------------------------------
    # 3. Micro-constants for the explicit two-compartment form (NONMEM
    #    ADVAN4: K = CL/V2, K23 = Q/V2, K32 = Q/V3).
    # ----------------------------------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ----------------------------------------------------------------------
    # 4. ODE system. Two-compartment disposition with a first-order depot
    #    arm and a parallel zero-order input straight into the central
    #    compartment.
    #
    #    DOSING. Each extended-release administration is TWO parallel dose
    #    records in the event table at the same time and the same nominal
    #    amount:
    #      - one to cmt = "depot"   (first-order arm, scaled by f(depot))
    #      - one to cmt = "central" with rate = -2 (modelled duration; the
    #        zero-order arm, scaled by f(central) over dur(central))
    #    An immediate-release administration needs only the depot record,
    #    because alpha and hence f(central) are zero.
    # ----------------------------------------------------------------------
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # ----------------------------------------------------------------------
    # 5. Bioavailability, zero-order duration and lag time.
    #    Control stream: F1 = TBIO*(1-INFRAC), F2 = TBIO*INFRAC, D2 = INDUR,
    #    ALAG1 = ALAG2 = LAG.
    # ----------------------------------------------------------------------
    f(depot)     <- ftot * (1 - alpha)
    f(central)   <- ftot * alpha
    dur(central) <- d2
    alag(depot)   <- tlag
    alag(central) <- tlag

    # ----------------------------------------------------------------------
    # 6. Observation. Doses are in mg and vc is in L, so central/vc is mg/L;
    #    the source reports plasma concentrations in ng/mL, hence the factor
    #    of 1000 (control stream: CP = A(2)/V2 on a dataset in consistent
    #    ng/mL units).
    # ----------------------------------------------------------------------
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
