Moein_2022_etrolizumab <- function() {
  description <- "Two-compartment population PK model for etrolizumab with first-order SC absorption and time-decreasing clearance in adults with moderately-to-severely active ulcerative colitis (Moein 2022)"
  reference <- "Moein A, Lu T, Jonsson S, et al. Population pharmacokinetic analysis of etrolizumab in patients with moderately-to-severely active ulcerative colitis. CPT Pharmacometrics Syst Pharmacol. 2022;11(9):1244-1255. doi:10.1002/psp4.12846"
  vignette <- "Moein_2022_etrolizumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric power scaling (WT / 70)^exponent on CL, Q, Vc, Vp (Table 4 notes a, b).",
      source_name        = "WT"
    ),
    ALB = list(
      description        = "Baseline serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Exponential effect on CL: exp(theta * (ALB - 41)). Reference 41 g/L is the median of the model-development cohort (Table 3) and the reference-patient value used in Figure 1.",
      source_name        = "ALB"
    ),
    CRP = list(
      description        = "Baseline C-reactive protein",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Exponential effect on CL: exp(theta * (CRP - 4.23)). Reference 4.23 mg/L is the reference-patient value reported in the Figure 1 forest plot caption (close to the Table 3 median of 4.31 mg/L). Standard CRP assay (typical of moderate-to-severe IBD populations where baseline CRP is well above the hs-CRP sensitivity range); the canonical general-scope CRP covariate covers both standard and high-sensitivity assays.",
      source_name        = "CRP"
    ),
    ADA_TITER = list(
      description        = "Time-varying antidrug antibody titer",
      units              = "(titer units)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Exponential effect on CL: exp(theta * ADA_TITER). Time-varying covariate -- use the actual titer at the time of each observation; paper imputes missing baseline as 0 (ADA negative) and postdose missingness by LOCF / NOCB. ADA-negative samples are encoded as ADA_TITER = 0 (American-spelling linear-titer convention, appropriate for exp(theta * ADA_TITER) effects); distinct from the reciprocal-dilution convention used in Jackson_2022_ixekizumab where ADA_TITER = 1 for negatives so log(1) = 0 cancels a log-linear effect.",
      source_name        = "ADAT"
    ),
    PRIOR_TNF = list(
      description        = "Prior anti-TNF inhibitor therapy indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no prior anti-TNF)",
      notes              = "Multiplicative effect on CL: (1 + theta * PRIOR_TNF).",
      source_name        = "PRIOR_TNF"
    ),
    DISEXT_EP = list(
      description        = "Disease extension: extensive colitis / pancolitis (vs. left-sided colitis)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (left-sided colitis; implies DISEXT_OTHER = 0 as well)",
      notes              = "Multiplicative effect on CL. Reference category left-sided colitis, the most common group (Table 3).",
      source_name        = "DISEXT"
    ),
    DISEXT_OTHER = list(
      description        = "Disease extension: other (neither left-sided colitis nor extensive/pancolitis)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (left-sided colitis)",
      notes              = "Multiplicative effect on CL. Mutually exclusive with DISEXT_EP. Only 2% of the development data fell in this group; the effect (+18%) carries large uncertainty (RSE 46.9%).",
      source_name        = "DISEXT"
    )
  )

  population <- list(
    n_subjects     = 1263,
    n_studies      = 5,
    age_range      = "18-79 years",
    age_median     = "38 years",
    weight_range   = "38.0-216 kg",
    weight_median  = "72.2 kg",
    sex_female_pct = 42,
    disease_state  = "Moderately-to-severely active ulcerative colitis; TNF-naive and TNF inadequate responders",
    dose_range     = "SC: 105 mg Q4W (phase III), 315-420 mg loading (phase II), 0.5-3 mg/kg Q4W (phase I); IV: 0.3-10 mg/kg single dose, 4 mg/kg Q4W (phase I)",
    regions        = "Multinational pooled phase I-III studies",
    ada_positive_pct = 23,
    prior_tnf_pct    = 45,
    disease_extension_pct = c(left_sided = 53, extensive_pancolitis = 42, other = 2, missing = 3),
    notes          = "Baseline demographics from Moein 2022 Table 3 (model-development cohort, N = 1263 with non-missing covariates). Contributing studies: ABS4262g (NCT00694980, phase I), EUCALYPTUS (NCT01336465, phase II), HIBISCUS I/II (NCT02163759/NCT02171429, phase III), HICKORY (NCT02100696, phase III), LAUREL (NCT02165215, phase III). GARDENIA (NCT02136069) was held out for external validation. Baseline medians: ALB 41 g/L, CRP 4.31 mg/L, fecal calprotectin 1500 ug/g."
  )

  ini({
    # Structural PK parameters -- reference values for a 70 kg adult at the first dose (TSFD = 0).
    # Source: Moein 2022 Table 4 (final population PK parameter estimates).
    lka      <- log(0.193);    label("First-order SC absorption rate constant (ka, 1/day)")                      # Table 4: ka = 0.193 /day
    lcl      <- log(0.260);    label("Baseline clearance at TSFD = 0 for 70 kg adult (CL, L/day)")               # Table 4: CL = 0.260 L/day
    lvc      <- log(2.61);     label("Central volume of distribution for 70 kg adult (Vc, L)")                   # Table 4: Vc = 2.61 L
    lvp      <- log(1.77);     label("Peripheral volume of distribution for 70 kg adult (Vp, L)")                # Table 4: Vp = 1.77 L
    lq       <- log(0.449);    label("Intercompartmental clearance for 70 kg adult (Q, L/day)")                  # Table 4: Q = 0.449 L/day
    logitfdepot  <- logit(0.712);  label("SC bioavailability (F, logit scale)")                                  # Table 4: F = 0.712

    # Time-dependent CL parameters (Equation 1, Moein 2022).
    logitmaxred  <- logit(0.263);  label("Maximum fractional reduction of CL over time (Maxred, logit scale)")   # Table 4: Maxred = 0.263
    lonset       <- log(4.81);     label("Half-life of the time-dependent CL change (Onset, weeks)")             # Table 4: Onset = 4.81 weeks

    # Allometric exponents on body weight (reference 70 kg; Table 4 notes a, b).
    allo_cl  <- 0.872;  label("Allometric exponent of WT on CL and Q")                                            # Table 4: WT on CL/Q = 0.872
    allo_v   <- 0.788;  label("Allometric exponent of WT on Vc and Vp")                                           # Table 4: WT on Vc/Vp = 0.788

    # Continuous covariate effects on CL (exponential form per Table 4 note d: CovEff = exp(theta * (Cov - Cov_ref))).
    e_alb_cl   <- -0.0314;   label("Albumin effect on CL (exp(theta * (ALB - 41)))")                             # Table 4: Albumin on CL = -0.0314
    e_crp_cl   <-  0.00458;  label("CRP effect on CL (exp(theta * (CRP - 4.23)))")                               # Table 4: CRP on CL = 0.00458
    e_adat_cl  <-  0.0365;   label("ADA titer effect on CL (exp(theta * ADA_TITER))")                            # Table 4: ADAT on CL = 0.0365

    # Categorical covariate effects on CL (multiplicative form per Table 4 note g: 1 + theta * indicator).
    e_priortnf_cl <- 0.0490; label("Prior anti-TNF fractional change in CL vs. no prior anti-TNF")                # Table 4: Prior TNF on CL = 0.0490
    e_extpan_cl   <- 0.0816; label("Extensive/pancolitis fractional change in CL vs. left-sided colitis")         # Table 4: Extensive/pancolitis on CL = 0.0816
    e_othext_cl   <- 0.181;  label("Other disease extension fractional change in CL vs. left-sided colitis")      # Table 4: Other disease extension on CL = 0.181

    # IIV (log-normal, variance on the estimation scale). omega^2 = log(1 + CV^2) from reported CVs.
    etalcl ~ log(1 + 0.243^2)  # Table 4: IIV CL CV = 0.243
    etalvc ~ log(1 + 0.252^2)  # Table 4: IIV Vc CV = 0.252
    etalvp ~ log(1 + 0.262^2)  # Table 4: IIV Vp CV = 0.262

    # IIV for F and Maxred is normally distributed on the logit-transformed scale.
    # Variance on the logit scale = (reported logit-SD)^2.
    etalogitfdepot  ~ 0.733^2   # Table 4: IIV F logit-SD = 0.733
    etalogitmaxred  ~ 0.597^2   # Table 4: IIV Maxred logit-SD = 0.597

    # Residual error (combined additive + proportional; Phase III typical values).
    propSd <- 0.196;  label("Proportional residual error (CV, fraction)")       # Table 4: proportional residual error CV = 0.196
    addSd  <- 0.427;  label("Additive residual error (ug/mL)")                  # Table 4: additive residual error SD = 0.427 ug/mL
  })

  model({
    # Time-dependent CL (Equation 1 of Moein 2022):
    #   CL_TSFD,j = CL_TSFD=0 * (1 - Maxred * (1 - exp(-log(2) / (Onset * 7) * (TSFD_j - TAD_j))))
    # TSFD_j - TAD_j is the administration time of the most-recent SC dose, relative to the first dose.
    # The model assumes SC dosing into depot; with tad(depot), t_since_first_dose = time in an rxode2 simulation
    # that starts at the first dose (time = 0), and time - tad(depot) gives the time of the most-recent dose.
    tdose        <- time - tad(depot)
    maxred       <- expit(logitmaxred + etalogitmaxred)
    onset        <- exp(lonset)
    td_cl_factor <- 1 - maxred * (1 - exp(-log(2) / (onset * 7) * tdose))

    # Covariate effects on CL (Table 4 notes d and g).
    alb_cl      <- exp(e_alb_cl  * (ALB - 41))
    crp_cl      <- exp(e_crp_cl  * (CRP - 4.23))
    adat_cl     <- exp(e_adat_cl * ADA_TITER)
    priortnf_cl <- 1 + e_priortnf_cl * PRIOR_TNF
    disext_cl   <- 1 + e_extpan_cl * DISEXT_EP + e_othext_cl * DISEXT_OTHER

    # Individual PK parameters with allometric weight scaling and covariate effects.
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (WT / 70)^allo_cl *
      td_cl_factor * alb_cl * crp_cl * adat_cl * priortnf_cl * disext_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^allo_v
    vp <- exp(lvp + etalvp) * (WT / 70)^allo_v
    q  <- exp(lq)           * (WT / 70)^allo_cl

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    f(depot) <- expit(logitfdepot + etalogitfdepot)

    # Concentration: dose in mg / volume in L = mg/L = ug/mL.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
