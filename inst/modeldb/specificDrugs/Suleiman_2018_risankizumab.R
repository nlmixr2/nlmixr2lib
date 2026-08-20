Suleiman_2018_risankizumab <- function() {
  description <- "Two-compartment population PK model of risankizumab (anti-IL-23 mAb) with first-order SC absorption in subjects with moderate-to-severe plaque psoriasis and moderate-to-severe Crohn's disease (Suleiman 2018 phase I-II)"
  reference <- "Suleiman AA, Khatri A, Minocha M, Othman AA. Population Pharmacokinetics of the Interleukin-23 Inhibitor Risankizumab in Subjects with Psoriasis and Crohn's Disease: Analyses of Phase I and II Trials. Clin Pharmacokinet. 2018;57(10):1259-1270. doi:10.1007/s40262-018-0704-z"
  vignette <- "Suleiman_2018_risankizumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Derived mechanically; verified = FALSE means it has
  # NOT been checked against the source paper.
  compartmentData <- list(
    depot       = list(analyte = "risankizumab", units = "mg", specimen = "administration site", verified = FALSE),
    central     = list(analyte = "risankizumab", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "risankizumab", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL, Q, Vc, and Vp; normalized as WT/70 per Suleiman 2018 Eq. 3 (reference 70 kg explicitly stated). One exponent (0.93) is shared between CL and Q; another (0.99) is shared between Vc and Vp.",
      source_name        = "WT"
    ),
    ALB = list(
      description        = "Baseline serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL; normalized as ALB/40 per Eq. 4 (reference 40 g/L, the overall cohort median from Table 2; verified by back-calculating the paper's typical-subject CL: 0.30 * (90/70)^0.93 * (42/40)^-1.54 = 0.351 L/day matches the paper's 0.35 L/day for a typical 90 kg psoriasis subject with 42 g/L albumin). Source uses SI units (g/L); convert g/dL to g/L by x10 if needed. Baseline-only, time-fixed per subject (paper evaluates ALB as a baseline covariate; ADAs were tested as time-varying, ALB was not).",
      source_name        = "ALB"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 272L,
    n_observations = 2831L,
    n_studies      = 3L,
    age_range      = "19-74 years (median 42.5)",
    age_median     = "42.5 years",
    weight_range   = "36-138 kg (median 79)",
    weight_median  = "79 kg",
    sex_female_pct = 44.1,
    race_ethnicity = c(White = 89.0, Hispanic = 7.0, Black = 2.6, Asian = 0.7, "American Indian/Alaska Native" = 0.7),
    disease_state  = "Pooled moderate-to-severe plaque psoriasis (n=157 from Studies 1 and 2, phase I single-ascending-dose and phase II multiple-dose) and moderate-to-severe Crohn's disease (n=115 from Study 3, phase II two-stage IV induction + SC maintenance). Baseline PASI median 16.2 (psoriasis); baseline CDAI median 297 (Crohn's).",
    dose_range     = "Study 1 (psoriasis, phase I): single dose 0.01-5 mg/kg IV or 0.25-1 mg/kg SC. Study 2 (psoriasis, phase II): 18 mg SC single dose, or 90 or 180 mg SC at weeks 0, 4, and 16. Study 3 (Crohn's, phase II): 200 or 600 mg IV every 4 weeks (induction) followed by 180 mg SC every 8 weeks (maintenance).",
    regions        = "Multi-center global (three studies).",
    notes          = "Baseline demographics from Suleiman 2018 Table 2 (n=272 pooled). 152/272 (55.9%) male, 120/272 (44.1%) female. Baseline albumin all-subject median 40 g/L (range 24-51; psoriasis median 42, Crohn's median 37). Baseline CRP all-subject median 4.3 mg/L. 30/272 (11%) developed treatment-emergent ADAs. Reference covariate values used in the covariate power terms: WT = 70 kg (explicitly stated in Sect. 2.3.2 / Eq. 3), ALB = 40 g/L (overall median, back-calculated against typical-subject CL values reported in the Abstract and Results Sect. 3.2). Below-quantification samples (169/3000, 5.63%) were excluded from the fit."
  )

  ini({
    # Structural parameters from Suleiman 2018 Table 3 (final phase I-II population PK model).
    # Typical values for the reference subject (70 kg, ALB 40 g/L).
    lcl     <- log(0.30);  label("Clearance (CL, L/day)")                              # Suleiman 2018 Table 3
    lvc     <- log(5.66);  label("Central volume of distribution (Vc, L)")             # Suleiman 2018 Table 3
    lq      <- log(0.33);  label("Intercompartmental clearance (Q, L/day)")            # Suleiman 2018 Table 3
    lvp     <- log(3.43);  label("Peripheral volume of distribution (Vp, L)")          # Suleiman 2018 Table 3
    lka     <- log(0.18);  label("First-order SC absorption rate (ka, 1/day)")         # Suleiman 2018 Table 3
    lfdepot <- log(0.72);  label("Absolute SC bioavailability (F, fraction)")          # Suleiman 2018 Table 3

    # Covariate effects (Suleiman 2018 Table 3; power exponents on covariates
    # normalized to their reference values, per Eq. 3 and Eq. 4). One exponent
    # is shared between CL and Q; another is shared between Vc and Vp.
    e_wt_cl_q  <-  0.93;   label("Shared power exponent of body weight on CL and Q (unitless)")   # Suleiman 2018 Table 3
    e_wt_vc_vp <-  0.99;   label("Shared power exponent of body weight on Vc and Vp (unitless)")  # Suleiman 2018 Table 3
    e_alb_cl   <- -1.54;   label("Power exponent of serum albumin on CL (unitless)")              # Suleiman 2018 Table 3

    # IIV. Table 3 reports %CV for the log-normal IIV entries with
    #   %CV = SQRT[exp(omega^2) - 1] * 100  (Table 3 footnote b),
    # so back-transform: omega^2 = log(CV^2 + 1).
    #   CL  37%  -> omega^2 = log(1 + 0.37^2) = 0.12832
    #   Vc  55%  -> omega^2 = log(1 + 0.55^2) = 0.26432
    #   Vp  35%  -> omega^2 = log(1 + 0.35^2) = 0.11556
    #   Ka  80%  -> omega^2 = log(1 + 0.80^2) = 0.49470
    # Suleiman 2018 Sect. 3.2 states "introducing a correlation between CL and Vc
    # reduced the OFV by 135 points" but Table 3 does NOT report the covariance
    # magnitude. The value is expected to be in the Electronic Supplementary
    # Material (ESM), which was not on disk during extraction.
    # Non-paper provenance: carried from sibling Suleiman_2019_risankizumab.R
    # (Suleiman AA is lead author on both; the 2018 phase I-II cohort is a
    # strict subset of the 2019 integrated phase I-III cohort). The 2019 model
    # reports a CL-Vc correlation of 39%; propagating that same correlation
    # onto the 2018 variances gives covariance = 0.39 * sqrt(0.12832 * 0.26432)
    # = 0.07182. See vignette Assumptions and deviations for the provenance
    # rationale and the operator sidecar exchange (request-001/002) that
    # authorised this substitution.
    etalcl + etalvc ~ c(0.12832,
                        0.07182, 0.26432)                                              # Suleiman 2018 Table 3 (IIV CL 37%, IIV Vc 55%); CL-Vc covariance carried from Suleiman 2019
    etalvp ~ 0.11556                                                                    # Suleiman 2018 Table 3 (IIV Vp 35%)
    etalka ~ 0.49470                                                                    # Suleiman 2018 Table 3 (IIV Ka 80%)

    # Residual error: proportional error model per Sect. 3.2. The magnitude is
    # NOT in Suleiman 2018 Table 3 (only IIV CVs and structural point estimates
    # are tabulated) and is expected to be in the unavailable ESM.
    # Non-paper provenance: carried from sibling Suleiman_2019_risankizumab.R,
    # which reports 19% CV. Same authorship group, same drug, 2018 phase I-II
    # cohort is a strict subset of the 2019 pooled cohort. See vignette
    # Assumptions and deviations.
    propSd <- 0.19;  label("Proportional residual error (SD, fraction); carried from Suleiman 2019")  # non-paper provenance: sibling Suleiman_2019_risankizumab
  })
  model({
    # Individual PK parameters. Reference subject: 70 kg, ALB 40 g/L. Covariate
    # forms per Suleiman 2018 Eq. 3 (WT on CL, Q, Vc, Vp with two shared
    # exponents) and Eq. 4 (ALB on CL).
    cl <- exp(lcl + etalcl) *
      (WT  / 70)^e_wt_cl_q *
      (ALB / 40)^e_alb_cl

    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc_vp
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vc_vp
    q  <- exp(lq)           * (WT / 70)^e_wt_cl_q
    ka <- exp(lka + etalka)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                k12 * central - k21 * peripheral1

    # Absolute SC bioavailability applied to the depot dose.
    f(depot) <- exp(lfdepot)

    # Concentration: dose in mg, volume in L -> mg/L = ug/mL.
    Cc <- central / vc

    Cc ~ prop(propSd)
  })
}
