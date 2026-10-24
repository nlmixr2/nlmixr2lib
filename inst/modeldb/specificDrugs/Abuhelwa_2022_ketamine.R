Abuhelwa_2022_ketamine <- function() {
  description <- "Joint parent-metabolite population PK model of racemic ketamine and norketamine after i.v., s.c. and i.m. ketamine in adults with treatment-refractory depression (Abuhelwa 2022). Ketamine has two-compartment disposition with first-order absorption shared by the s.c. and i.m. routes (common ka and common absolute bioavailability of 64.4 percent); all ketamine clearance forms norketamine (fraction metabolised assumed 1, so the norketamine clearances and volumes are apparent values divided by fm), which has its own two-compartment disposition. All clearances and volumes are scaled allometrically to 70 kg (exponent 0.75 for clearances, 1 for volumes). Residual error is proportional for ketamine (additive term fixed near zero) and combined additive plus proportional for norketamine."
  reference <- "Abuhelwa AY, Somogyi AA, Loo CK, Glue P, Barratt DT, Foster DJR. Population Pharmacokinetics and Pharmacodynamics of the Therapeutic and Adverse Effects of Ketamine in Patients With Treatment-Refractory Depression. Clin Pharmacol Ther. 2022;112(3):720-729. doi:10.1002/cpt.2640."
  vignette <- "Abuhelwa_2022_ketamine"
  units <- list(time = "h", dosing = "ug", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric scaling referenced to 70 kg on every clearance (exponent 0.75: CL, Q, CL/fm, Q/fm) and every volume (exponent 1: V1, V2, V1/fm, V2/fm) of ketamine and norketamine (Methods 'General modeling approach'; supplementary NONMEM $PK, WTS = 70).",
      source_name = "WEIGHT"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "ketamine", units = "ug", specimen = "administration site", verified = TRUE),
    central = list(analyte = "ketamine", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "ketamine", units = "ug", specimen = "plasma", verified = TRUE),
    central_nk = list(
      analyte = "norketamine",
      units = "ug (ketamine mass equivalents)",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1_nk = list(
      analyte = "norketamine",
      units = "ug (ketamine mass equivalents)",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 21L,
    n_studies = 2L,
    age_range = "29-69 years",
    age_median = "50 years (mean 49.5)",
    weight_range = "not reported",
    sex_female_pct = 33.3,
    disease_state = "Treatment-refractory major depression (MADRS >= 20 before each treatment)",
    dose_range = "Racemic ketamine 0.1-0.5 mg/kg single doses in a within-subject dose-titration design, by i.v., s.c. or i.m. injection; doses at least 1 week apart",
    regions = "Australia (University of New South Wales, Sydney)",
    genotype = "CYP2B6 *6/*6 : *1/*6 or *5 : *1/*1 = 2 : 8 : 11",
    n_observations = "278 ketamine and 278 norketamine plasma concentrations (Table S1); 14 percent of concentrations were below the 0.25 ng/mL LLOQ and were excluded (M1 method)",
    notes = "Demographics from Abuhelwa 2022 Table 1 (7 female, 14 male). Data are from 61 ketamine treatment sessions of the Loo 2016 multiple-route dose-titration crossover study, including the Lai 2014 i.v. pilot. Body weight was used as a covariate but its distribution is not reported."
  )

  ini({
    # Ketamine (Table 2; supplementary NONMEM $THETA). Typical values for 70 kg.
    lka <- log(6.53); label("First-order absorption rate constant, common to s.c. and i.m. (1/h)") # Table 2 'KA ketamine' = 6.53 h-1 (RSE 19.1%); one THETA for both routes
    lfdepot <- log(0.644); label("Absolute bioavailability of s.c. and i.m. ketamine (fraction)") # Table 2 'F extravascular' = 0.644 (RSE 11.5%); F2 = F3 in the supplementary code
    lcl <- log(69.6); label("Ketamine clearance at 70 kg (L/h)") # Table 2 'CL ketamine' = 69.6 L/h (RSE 10.5%)
    lvc <- log(79.3); label("Ketamine central volume of distribution at 70 kg (L)") # Table 2 'V1 ketamine' = 79.3 L (RSE 18.7%)
    lq <- log(121); label("Ketamine intercompartmental clearance at 70 kg (L/h)") # Table 2 'Q ketamine' = 121 L/h (RSE 14.3%)
    lvp <- log(87.4); label("Ketamine peripheral volume of distribution at 70 kg (L)") # Table 2 'V2 ketamine' = 87.4 L (RSE 14.8%)

    # Norketamine apparent parameters (divided by fm, which is assumed to be 1;
    # Methods 'Pharmacokinetic modeling').
    lcl_nk <- log(30.3); label("Norketamine apparent clearance CL/fm at 70 kg (L/h)") # Table 2 'CL norketamine/fm' = 30.3 L/h (RSE 31.2%)
    lvc_nk <- log(73.6); label("Norketamine apparent central volume V1/fm at 70 kg (L)") # Table 2 'V1 norketamine/fm' = 73.6 L (RSE 12.4%)
    lq_nk <- log(54.2); label("Norketamine apparent intercompartmental clearance Q/fm at 70 kg (L/h)") # Table 2 'Q norketamine/fm' = 54.2 L/h (RSE 11.3%)
    lvp_nk <- log(339); label("Norketamine apparent peripheral volume V2/fm at 70 kg (L)") # Table 2 'V2 norketamine/fm' = 339 L (RSE 34.2%)

    # Allometric exponents, fixed (Methods 'General modeling approach').
    e_wt_cl_q <- fixed(0.75); label("Allometric exponent of WT on ketamine CL and Q (unitless)") # Methods: 'clearances scaled to the exponent of 0.75'
    e_wt_vc_vp <- fixed(1); label("Allometric exponent of WT on ketamine V1 and V2 (unitless)") # Methods: 'volume to 1'
    e_wt_cl_nk <- fixed(0.75); label("Allometric exponent of WT on norketamine CL/fm (unitless)") # Methods: 'clearances scaled to the exponent of 0.75'
    e_wt_q_nk <- fixed(0.75); label("Allometric exponent of WT on norketamine Q/fm (unitless)") # Methods: 'clearances scaled to the exponent of 0.75'
    e_wt_vc_nk <- fixed(1); label("Allometric exponent of WT on norketamine V1/fm (unitless)") # Methods: 'volume to 1'
    e_wt_vp_nk <- fixed(1); label("Allometric exponent of WT on norketamine V2/fm (unitless)") # Methods: 'volume to 1'

    # IIV: Table 2 BSV %CV = sqrt(omega2) x 100, confirmed by the supplementary
    # $OMEGA values (e.g. sqrt(0.0341) = 0.185).
    etalcl ~ 0.0341 # Table 2 BSV 18.5%CV; supplementary $OMEGA PPV_CLK = 0.0341
    etalvc ~ 0.219 # Table 2 BSV 46.8%CV; supplementary $OMEGA PPV_V1K = 0.219
    etalvp ~ 0.182 # Table 2 BSV 42.7%CV; supplementary $OMEGA PPV_V2K = 0.182
    etalcl_nk ~ 0.424 # Table 2 BSV 65.1%CV; supplementary $OMEGA PPV_CLNK = 0.424
    etalvc_nk ~ 0.217 # Table 2 BSV 46.6%CV; supplementary $OMEGA PPV_V1NK = 0.217

    # Residual error: Y = IPRED + sqrt((IPRED * prop)^2 + add^2) * EPS, SIGMA 1 FIX.
    propSd <- 0.294; label("Proportional residual error, ketamine (fraction)") # Table 2 'RUVPROP-ketamine' = 0.294 (RSE 9.6%)
    addSd <- fixed(0.00001); label("Additive residual error, ketamine (ng/mL)") # Table 2 'RUVADD-ketamine' = 0.00001 Fix
    propSd_nk <- 0.257; label("Proportional residual error, norketamine (fraction)") # Table 2 'RUVPROP-norketamine' = 0.257 (RSE 18.6%)
    addSd_nk <- 3.88; label("Additive residual error, norketamine (ng/mL)") # Table 2 'RUVADD-norketamine' = 3.88 ng/mL (RSE 14.7%)
  })

  model({
    # Individual parameters with allometric scaling to 70 kg (supplementary $PK).
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl_q
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc_vp
    q <- exp(lq) * (WT / 70)^e_wt_cl_q
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vc_vp

    cl_nk <- exp(lcl_nk + etalcl_nk) * (WT / 70)^e_wt_cl_nk
    vc_nk <- exp(lvc_nk + etalvc_nk) * (WT / 70)^e_wt_vc_nk
    q_nk <- exp(lq_nk) * (WT / 70)^e_wt_q_nk
    vp_nk <- exp(lvp_nk) * (WT / 70)^e_wt_vp_nk

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    kel_nk <- cl_nk / vc_nk
    k12_nk <- q_nk / vc_nk
    k21_nk <- q_nk / vp_nk

    # i.v. doses go to central; s.c. and i.m. doses go to depot (the NONMEM
    # code keeps separate s.c. and i.m. depots but gives them the same KA and
    # F, so one depot is equivalent). All ketamine clearance forms
    # norketamine in ketamine mass units (fm = 1, no molecular-weight
    # conversion), exactly as in the supplementary $DES.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(central_nk) <- kel * central - kel_nk * central_nk - k12_nk * central_nk + k21_nk * peripheral1_nk
    d/dt(peripheral1_nk) <- k12_nk * central_nk - k21_nk * peripheral1_nk

    f(depot) <- exp(lfdepot)

    # ug / L = ng/mL
    Cc <- central / vc
    Cc_nk <- central_nk / vc_nk

    Cc ~ add(addSd) + prop(propSd)
    Cc_nk ~ add(addSd_nk) + prop(propSd_nk)
  })
}
