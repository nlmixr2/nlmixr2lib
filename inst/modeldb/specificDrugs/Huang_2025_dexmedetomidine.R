Huang_2025_dexmedetomidine <- function() {
  description <- paste(
    "Two-compartment population PK model for dexmedetomidine nasal spray in",
    "Chinese adults, pooled across one phase I study in healthy volunteers",
    "(intravenous and intranasal) and one phase III study in patients",
    "undergoing elective abdominal surgery (Huang 2025; n = 196 subjects,",
    "1,225 plasma concentrations). First-order absorption from a nasal depot",
    "with an absorption lag time and estimated bioavailability (65.3%), and",
    "linear elimination. Body weight enters as theory-based allometric",
    "scaling with exponents fixed a priori at 0.75 on CL and Q and 1 on Vc",
    "and Vp, normalised to the 60 kg cohort median. The single retained",
    "covariate beyond size is health status: the absorption rate constant in",
    "healthy volunteers is 2.05-fold that in patients (equivalently, patient",
    "KA is 49% of healthy-volunteer KA), which the authors attribute partly",
    "to the sparse phase III sampling during the absorption phase.",
    "Inter-individual variability is estimated on CL and KA only. Residual",
    "error is combined proportional plus additive. The paper's companion",
    "exposure-response analysis is a logistic regression of the probability",
    "of a Ramsay Sedation Scale score >= 3 within 45 min on Cmax; that is a",
    "non-ODE statistical regression on a derived exposure metric, so",
    "following the Yin_2021_pexidartinib.R precedent it is reproduced in the",
    "validation vignette rather than encoded here."
  )
  reference <- paste(
    "Huang Y, Xu S, Djebli N, Jiang H, Yang G (2025). Population",
    "pharmacokinetic modeling of dexmedetomidine nasal spray in Chinese",
    "adults. Frontiers in Pharmacology 16:1662364.",
    "doi:10.3389/fphar.2025.1662364.",
    sep = " "
  )
  vignette <- "Huang_2025_dexmedetomidine"
  units    <- list(time = "h", dosing = "ug", concentration = "pg/mL")

  compartmentData <- list(
    depot = list(
      analyte = "dexmedetomidine", units = "ug",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "dexmedetomidine", units = "ug",
      specimen = "plasma", verified = TRUE
    ),
    peripheral1 = list(
      analyte = "dexmedetomidine", units = "ug",
      specimen = "tissue", verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Drives theory-based allometric scaling of CL and Q (exponent 0.75)",
        "and of Vc and Vp (exponent 1), both fixed a priori rather than",
        "estimated. Huang 2025 Results 'Pharamcokinetic modeling': 'with BW",
        "incorporated via theory-based allometric scaling on distribution and",
        "elimination (i.e., on CL, Vc, Vp, and Q)'; Discussion: 'with",
        "exponents of 0.75 and 1 for clearances and volumes, respectively'.",
        "REFERENCE WEIGHT = 60 kg. The paper never writes the normalising",
        "weight into a printed equation, but it is sourced rather than",
        "assumed: Methods Eq. 6 defines the power model as",
        "Pi = PTV * (COV / COVmedian)^theta with 'COVmedian is the median or",
        "typical value of the covariate', and Table 3 gives the total-cohort",
        "median BW = 60 kg (Phase I 57.65, Phase III 60). The competing",
        "reading is the 70 kg standard adult implied by the phrase",
        "'theory-based allometric scaling'; the two differ by 17% in every",
        "predicted concentration. Scored against the paper's own simulation",
        "outputs, 60 kg reproduces all eight Figure 4 points and the three",
        "pediatric mean Cmax values as arithmetic means to within a few",
        "percent. Those mean-scale anchors are by themselves DEGENERATE: 70 kg",
        "scored against simulated medians fits them equally well, because the",
        "mean/median ratio of the simulated Cmax distribution is itself close",
        "to 70/60. The tie is broken by the six pediatric target-attainment",
        "percentages printed in Results ('Pediatric exposure simulation'),",
        "which are insensitive to that degeneracy and select 60 kg: at 70 kg",
        "the 180 pg/mL attainment is biased high at all three strata, by",
        "about 6 percentage points on average, whereas at 60 kg the same",
        "bias is near zero. Figure 5 is the one published quantity",
        "pointing the other way, and its caption does not describe its own",
        "content -- see the validation vignette Errata for the full scoring",
        "tables and the Figure 5 defect."
      ),
      source_name        = "BW"
    ),
    DIS_HEALTHY = list(
      description        = "Healthy-participant cohort indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (patient undergoing elective abdominal surgery, phase III NCT04383418)",
      notes              = paste(
        "Source column 'state (HV vs. patients)'. Applied to the absorption",
        "rate constant only, in the Methods Eq. 7 piecewise categorical form",
        "Pi = PTV if COV = type1, PTV * (1 + theta) if COV = type2. Table 4",
        "gives 'state on KA' = 1.05, so KA_healthy = 0.523 * (1 + 1.05) =",
        "1.072 /h and KA_patient = 0.523 /h. The orientation (patients are",
        "the reference, type1) is fixed by two independent statements in the",
        "paper: the Abstract, 'The absorption rate in the patients from phase",
        "III study was approximately 49% of that in the HV from phase I",
        "study' (0.523 / 1.072 = 48.8%), and the Discussion, 'The KA in",
        "healthy volunteers was approximately 1 h-1'. A multiplicative",
        "reading (KA * 1.05) is incompatible with both. Independently",
        "confirmed against Figure 4, which plots simulated Cmax for both",
        "groups at 40/60/80/100 kg: the published healthy-to-patient Cmax",
        "ratio is 1.48-1.51 across all four weights, matching the 2.05-fold",
        "KA reading (simulated 1.44-1.54) and excluding the 1.05-fold reading",
        "(which would give ~1.02).",
        "The phase I cohort is entirely healthy and the phase III cohort",
        "entirely patients, so in this analysis the indicator is fully",
        "confounded with study, sampling density (rich vs sparse), and",
        "surgical/anaesthetic context; the authors note the absorption-phase",
        "sampling sparsity as a likely partial explanation and suggest the",
        "true difference 'might be less pronounced as compared to the",
        "estimated value'."
      ),
      source_name        = "state"
    )
  )

  # Covariates screened in the stepwise forward-inclusion / backward-
  # elimination procedure but NOT retained in the final model. Huang 2025
  # Methods: 'The covariates under investigation include state (HV vs.
  # patients), sex, body mass index (BMI), body weight (BW), body surface
  # area (BSA), and age.' Results: BSA and BMI were dropped a priori for
  # multicollinearity with BW; age, sex, and state were tested on KA and CL,
  # and only the state-on-KA effect survived. No point estimates are reported
  # for the rejected effects, so they are documented here rather than encoded.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Screened on KA and CL and not retained. Table 3 median 38 years",
        "(range 18-65); the phase III cohort is materially older (median 44)",
        "than the phase I cohort (median 22)."
      ),
      source_name = "AGE"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened on KA and CL and not retained. Table 3 reports 132 female",
        "and 64 male of 196 (67.3% female), with a higher female proportion",
        "in phase III (108/148) than phase I (24/48). Source column GENDER."
      ),
      source_name = "GENDER"
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = paste(
        "Excluded before stepwise screening because of multicollinearity",
        "with body weight (Results: 'Due to the multicollinearity between",
        "BW, BSA and BMI, BSA and BMI were not considered during the",
        "covariates' screening after introducing BW'). Table 3 median",
        "22.9 kg/m^2 (range 18.6-29.9)."
      ),
      source_name = "BMI"
    ),
    BSA = list(
      description = "Body surface area",
      units       = "m^2",
      type        = "continuous",
      notes       = paste(
        "Excluded before stepwise screening because of multicollinearity",
        "with body weight (same Results sentence as BMI). Table 3 median",
        "1.6 m^2 (range 1.32-2.19)."
      ),
      source_name = "BSA"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 196,
    n_studies      = 2,
    n_observations = 1225,
    age_range      = "18-65 years",
    age_median     = "38 years",
    weight_range   = "45.2-97 kg",
    weight_median  = "60 kg",
    sex_female_pct = 67.3,
    race_ethnicity = "Chinese (all subjects enrolled in China; the paper reports no race/ethnicity breakdown)",
    disease_state  = paste(
      "Healthy volunteers (phase I, n = 48) and adults undergoing elective",
      "abdominal surgery excluding liver surgery, requiring general",
      "anaesthesia, endotracheal intubation and mechanical ventilation",
      "(phase III, n = 148)"
    ),
    dose_range     = paste(
      "Phase I: 20 and 40 ug intravenous over 15 min; 150 ug and 100/20 ug",
      "intranasal. Phase III: 75 and 100 ug intranasal (single dose)"
    ),
    regions        = "China",
    notes          = paste(
      "Baseline demographics from Huang 2025 Table 3 (medians with min-max).",
      "Study designs from Table 1: phase I registrations CTR20191868 and",
      "CTR20171118 (three parts, intensive sampling to 10-24 h) and phase III",
      "NCT04383418 (sparse sampling, two samples per subject). The phase I",
      "study is reported in Kuang et al. 2022. Approximately 1% of samples",
      "were below the limit of quantification and were handled by the M1",
      "method (discarded)."
    )
  )

  ini({
    # ---- Structural parameters -------------------------------------------
    # Huang 2025 Table 4 'Final model / Estimates (RSE %)'. These are the
    # typical values at the 60 kg reference weight; KA is the PATIENT value
    # (DIS_HEALTHY = 0), see the covariateData note.
    lcl <- log(35.3);    label("Clearance CL (L/h) at WT = 60 kg")                                  # Table 4 'CL, L/h' = 35.3 (RSE 3.3%; bootstrap median 35.34, 95% CI 32.787-37.689)
    lvc <- log(21.5);    label("Central volume of distribution Vc (L) at WT = 60 kg")               # Table 4 'Vc, L' = 21.5 (RSE 6.3%; bootstrap median 21.51, 95% CI 19.128-25.007)
    lq  <- log(116);     label("Inter-compartmental clearance Q (L/h) at WT = 60 kg")               # Table 4 'Q, L/h' = 116 (RSE 5.8%; bootstrap median 115.58, 95% CI 103.304-129.962)
    lvp <- log(86.5);    label("Peripheral volume of distribution Vp (L) at WT = 60 kg")            # Table 4 'Vp, L' = 86.5 (RSE 3.01%; bootstrap median 86.42, 95% CI 80.384-92.26)
    lka <- log(0.523);   label("Absorption rate constant KA (1/h), patient reference cohort")       # Table 4 'KA, h -1' = 0.523 (RSE 9%; bootstrap median 0.53, 95% CI 0.437-0.629)

    ltlag   <- log(0.0592); label("Absorption lag time ALAG (h)")                                   # Table 4 'ALAG,h' = 0.0592 (RSE 8.6%; bootstrap median 0.0594, 95% CI 0.047-0.068)
    lfdepot <- log(0.653);  label("Bioavailability F1 of the nasal spray (fraction)")               # Table 4 'F1, %' = 0.653 (RSE 3.9%; bootstrap median 0.655, 95% CI 0.592-0.706)

    # ---- Allometric exponents, fixed a priori ----------------------------
    # Discussion: 'Body weight was first included as a covariate affecting the
    # volumes of distribution and clearances in the base model via theory-based
    # allometric scaling (i.e., with exponents of 0.75 and 1 for clearances and
    # volumes, respectively).' Not estimated -- no Table 4 row, no RSE.
    e_wt_cl_q  <- fixed(0.75); label("Allometric exponent of body weight on CL and Q (unitless)")   # Results 'Pharamcokinetic modeling' + Discussion paragraph 4 (WBE theory, West 1997/1999)
    e_wt_vc_vp <- fixed(1.0);  label("Allometric exponent of body weight on Vc and Vp (unitless)")  # Results 'Pharamcokinetic modeling' + Discussion paragraph 4 (Holford & Anderson 2017)

    # ---- Covariate effect ------------------------------------------------
    # Eq. 7 piecewise categorical, applied as KA * (1 + theta) for healthy
    # volunteers. See the DIS_HEALTHY covariateData note for the four-way
    # confirmation of the orientation and the additive-(1+theta) form.
    e_dis_healthy_ka <- 1.05; label("Healthy-volunteer effect on KA (Eq. 7 additive factor; KA_HV = KA * (1 + 1.05))")  # Table 4 'state on KA' = 1.05 (RSE 27.5%; bootstrap median 1.06, 95% CI 0.525-1.719)

    # ---- Inter-individual variability ------------------------------------
    # Eq. 1 is Pi = PTV * exp(eta_i), so omega is the SD of eta on the log
    # scale and the Table 4 'omega (CL/KA), %' rows are that SD as a percent.
    # Two checks confirm the SD reading over a variance reading:
    #  (a) a variance RSE cannot fall below sqrt(2/N) = sqrt(2/196) = 10.1%,
    #      but Table 4 reports RSE 9.5% for omega(KA);
    #  (b) both 95% CIs are exactly symmetric about the estimate
    #      (22.4 +/- 4.61, 92.3 +/- 17.21), i.e. the tabulated number is the
    #      estimated quantity itself rather than a nonlinear transform of it.
    # The supplement's run-acceptance criteria corroborate the scale by
    # discussing 'the standard error of Etas estimates' as a percent of the
    # estimate. IIV was estimated on CL and KA only.
    etalcl ~ 0.224^2   # Table 4 'omega (CL), %' = 22.4 (RSE 10.5%; 95% CI 17.79-27.01) -> variance 0.050176
    etalka ~ 0.923^2   # Table 4 'omega (KA), %' = 92.3 (RSE  9.5%; 95% CI 75.09-109.51) -> variance 0.851929

    # ---- Residual error --------------------------------------------------
    # Eq. 4 combined model: Yobs = Ypred * (1 + eps1) + eps2. Table 4 reports
    # the proportional term as a percent and the additive term in pg/mL, the
    # units the assay and the whole paper report concentrations in.
    propSd <- 0.277; label("Proportional residual error (fraction)")  # Table 4 'sigma (Prop), %'   = 27.7 (RSE  5.9%; 95% CI 24.51-30.89)
    addSd  <- 5.01;  label("Additive residual error (pg/mL)")         # Table 4 'sigma (Add), pg/mL' = 5.01 (RSE 16.7%; 95% CI 3.37-6.65)
  })

  model({
    # 1. Body-weight allometry. Eq. 6 power form normalised to COVmedian; the
    #    cohort median body weight is 60 kg (Table 3, Total column).
    allom_cl <- (WT / 60)^e_wt_cl_q
    allom_v  <- (WT / 60)^e_wt_vc_vp

    # 2. Individual parameters. IIV on CL and KA only (Results: 'The
    #    inter-individual variability (IIV) terms were estimated on CL and
    #    KA'). The health-status effect enters KA in the Eq. 7 additive form.
    cl <- exp(lcl + etalcl) * allom_cl
    vc <- exp(lvc)          * allom_v
    q  <- exp(lq)           * allom_cl
    vp <- exp(lvp)          * allom_v

    ka     <- exp(lka + etalka) * (1 + e_dis_healthy_ka * DIS_HEALTHY)
    tlag   <- exp(ltlag)
    fdepot <- exp(lfdepot)

    # 3. Micro-constants
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # 4. Two-compartment disposition with first-order absorption from the
    #    nasal depot. Intravenous doses (phase I parts 1) are given directly
    #    into `central`, where neither F1 nor the lag time applies; intranasal
    #    doses are given into `depot`.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # 5. Bioavailability and lag time apply to the nasal route only.
    f(depot)    <- fdepot
    alag(depot) <- tlag

    # 6. Plasma dexmedetomidine concentration. Dose is in ug and volumes are
    #    in L, so central/vc is ug/L = ng/mL; the paper reports concentrations
    #    in pg/mL throughout (Table 4 additive sigma, the 100 and 180 pg/mL
    #    exposure-response targets, Figures 3-5), and 1 ng/mL = 1000 pg/mL.
    Cc <- central / vc * 1000
    Cc ~ add(addSd) + prop(propSd)
  })
}
