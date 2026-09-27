Okamoto_2021_vedolizumab <- function() {
  description <- "Two-compartment population PK model for vedolizumab (humanised anti-alpha4-beta7 integrin IgG1 monoclonal antibody) with parallel linear and Michaelis-Menten elimination, a time-varying anti-vedolizumab-antibody titer effect and an Asian-race effect on linear clearance, in Asian and non-Asian adults with moderately-to-severely active ulcerative colitis or Crohn's disease (Okamoto 2021)."
  reference <- "Okamoto H, Dirks NL, Rosario M, Hori T, Hibi T. Population pharmacokinetics of vedolizumab in Asian and non-Asian patients with ulcerative colitis and Crohn's disease. Intest Res. 2021;19(1):95-105. doi:10.5217/ir.2019.09167 (PMC7873400)."
  vignette <- "Okamoto_2021_vedolizumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight (time-varying)",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying, carried forward from the last observation (Okamoto 2021 Methods 'Data Assembly'). Power-form effects on CLL (exponent 0.472, estimated) and Vc (0.466, estimated); fixed allometric exponents on Vp (1), Q (0.75) and Vmax (0.75). Reference 70 kg (Table 2 footnote).",
      source_name = "WT"
    ),
    ALB = list(
      description = "Serum albumin (time-varying)",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying, carried forward from the last observation. Power-form effect on CLL, (ALB / 40 g/L)^(-1.19). The paper reports albumin in g/dL with a reference of 4.0 g/dL; the canonical ALB column is in g/L, so the reference is written as 40 g/L (4.0 g/dL x 10).",
      source_name = "ALB"
    ),
    ADA_TITER = list(
      description = "Anti-vedolizumab antibody (AVA) titer by electrochemiluminescence assay (time-varying). Linear-titer zero-encoding convention: ADA_TITER = 0 for AVA-negative records (any titer below 10).",
      units = "(reciprocal dilution titer)",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying; the authors linearly interpolated titers between observed samples. AVA-positive is defined as a titer >= 10 and AVA-negative as a titer < 10 (Okamoto 2021 Results 'Final PK Model'). A record with ADA_TITER >= 10 switches the typical CLL from the AVA-negative value (0.165 L/day) to the AVA-positive value (0.246 L/day at a titer of 250) and applies the titer power term (ADA_TITER / 250)^0.0713; a record with ADA_TITER < 10 (including 0) takes the AVA-negative CLL with no titer term. Total binding antibodies, not neutralizing antibodies. Observed positive titers were typically 10-6,250.",
      source_name = "AVA titer"
    ),
    RACE_ASIAN = list(
      description = "Race indicator (1 = Asian, 0 = non-Asian)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-Asian)",
      notes = "Multiplicative effect on CLL of the form 1.10^RACE_ASIAN. Non-Asian pools White, Black, American Indian or Alaska Native, Native Hawaiian or other Pacific Islander, and Other (Okamoto 2021 Results 'Pharmacokinetic Analysis Dataset'). Baseline (not time-varying).",
      source_name = "Race (Asian vs non-Asian)"
    ),
    IBD_CD = list(
      description = "IBD diagnosis indicator (1 = Crohn's disease, 0 = ulcerative colitis)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (ulcerative colitis)",
      notes = "Multiplicative effect on CLL of the form 0.990^IBD_CD. Reference diagnosis UC (Table 2 footnote). Baseline (not time-varying).",
      source_name = "Diagnosis (UC vs CD)"
    ),
    OCC = list(
      description = "Integer-valued dosing-occasion indicator for inter-occasion variability on CLL",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = "Okamoto 2021 Methods 'Data Assembly' defines an occasion as a dosing interval with at least one associated vedolizumab PK measurement, but does not state the number of occasions. Ten occasions are encoded here, which covers the densest per-patient schedule in the pooled studies (the GEMINI 1/2 every-4-weeks maintenance arm, whose PK samples fall in 10 distinct dosing intervals; Supplementary Table 1). Decomposed inside model() into binary indicators oc1..oc10 multiplexing ten IOV etas on log-CLL that share one variance. Records with OCC outside 1..10 carry no IOV.",
      source_name = "OCC"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Present in the analysis dataset (Okamoto 2021 Table 1: median 36, range 17-79 years) but not included in the full covariate model."
    ),
    SEXF = list(
      description = "Sex (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      notes = "Present in the analysis dataset (47.3% female; Table 1) but not included in the full covariate model."
    )
  )

  compartmentData <- list(
    central = list(analyte = "vedolizumab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "vedolizumab", units = "mg", specimen = "serum", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 1933L,
    n_studies = 5L,
    n_observations = 16598L,
    age_range = "17-79 years",
    age_median = "36 years",
    weight_range = "28-172 kg (baseline)",
    weight_median = "67 kg (baseline)",
    albumin_range = "1.40-5.30 g/dL (baseline; median 3.80 g/dL)",
    sex_female_pct = 47.3,
    race_ethnicity = c(
      White = 76.9,
      Asian = 20.5,
      Black = 1.6,
      `American Indian or Alaska Native` = 0.3,
      `Native Hawaiian or other Pacific Islander` = 0.2,
      Other = 0.6
    ),
    disease_state = "Moderately-to-severely active ulcerative colitis (46.8%) or Crohn's disease (53.2%).",
    dose_range = "300 mg IV at weeks 0, 2 and 6, then 300 mg IV every 4 or every 8 weeks (phase 3); 150 or 300 mg IV on days 1, 15 and 43 (phase 1 CPH-001).",
    regions = "Japan (CPH-001, CCT-101, CCT-001; n = 224 Japanese patients); North America, Europe, Asia/Australia/Africa (GEMINI 1 and 2).",
    ada_positive_pct = 8.2,
    asian_n = 396L,
    notes = "Pooled phase 1 CPH-001 (Japanese UC, n = 9), phase 3 CCT-101 (Japanese UC, n = 152), CCT-001 (Japanese CD, n = 63), GEMINI 1 (UC, n = 743) and GEMINI 2 (CD, n = 966). Demographics from Okamoto 2021 Table 1 and Results 'Pharmacokinetic Analysis Dataset'. Placebo-only patients were excluded."
  )

  ini({
    # Structural parameters -- Okamoto 2021 Table 2 (medians of the Bayesian
    # posterior; NONMEM 7.3 MCMC with $PRIOR). Reference patient: 70 kg,
    # albumin 4 g/dL, non-Asian, UC, AVA-negative (titer < 10); AVA-positive
    # reference titer 250.
    lcl <- log(0.165); label("Linear clearance CLL for the reference AVA-negative patient (L/day)") # Table 2: AVA- CLL = 0.165 L/day (0.160, 0.170)
    lcl_adapos <- log(0.246); label("Linear clearance CLL for the reference AVA-positive patient at titer 250 (L/day)") # Table 2: AVA+ CLL = 0.246 L/day (0.222, 0.273)
    lvc <- log(3.16); label("Central volume of distribution Vc (L)") # Table 2: Vc = 3.16 L (3.11, 3.22)
    lvp <- log(1.84); label("Peripheral volume of distribution Vp (L)") # Table 2: Vp = 1.84 L (1.71, 1.99)
    lq <- log(0.161); label("Intercompartmental clearance Q (L/day)") # Table 2: Q = 0.161 L/day (0.150, 0.173)
    lvmax <- log(0.238); label("Maximum elimination rate of the Michaelis-Menten pathway (mg/day)") # Table 2: Vmax = 0.238 mg/day (0.191, 0.296)
    lkm <- log(0.851); label("Michaelis-Menten constant Km (ug/mL)") # Table 2: Km = 0.851 ug/mL (0.641, 1.150)

    # Covariate effects on CLL and Vc (Table 2).
    e_wt_cl <- 0.472; label("Weight exponent on CLL (unitless; reference 70 kg)") # Table 2: CLL ~ WT = 0.472 (0.400, 0.533)
    e_alb_cl <- -1.19; label("Albumin exponent on CLL (unitless; reference 4 g/dL = 40 g/L)") # Table 2: CLL ~ albumin = -1.19 (-1.27, -1.11)
    e_ada_titer_cl <- 0.0713; label("AVA-titer power exponent on CLL for AVA-positive records, (ADA_TITER/250)^e (unitless)") # Table 2: CLL ~ AVA+ = 0.0713 (0.0404, 0.1020)
    e_race_asian_cl <- 1.10; label("Asian-race multiplier on CLL, CLL * e^RACE_ASIAN (unitless)") # Table 2: CLL ~ race: Asian = 1.10 (1.06, 1.14)
    e_ibd_cd_cl <- 0.990; label("Crohn's-disease multiplier on CLL, CLL * e^IBD_CD (unitless)") # Table 2: CLL ~ diagnosis: CD = 0.990 (0.960, 1.020)
    e_wt_vc <- 0.466; label("Weight exponent on Vc (unitless; reference 70 kg)") # Table 2: Vc ~ WT = 0.466 (0.424, 0.509)

    # Fixed allometric exponents (Table 2 'Fixed').
    e_wt_vp <- fixed(1); label("Allometric exponent of WT on Vp (unitless)") # Table 2: Vp ~ WT = 1.00 Fixed
    e_wt_vmax <- fixed(0.75); label("Allometric exponent of WT on Vmax (unitless)") # Table 2: Vmax ~ WT = 0.750 Fixed
    e_wt_q <- fixed(0.75); label("Allometric exponent of WT on Q (unitless)") # Table 2: Q ~ WT = 0.750 Fixed

    # Interindividual variability -- Table 2, full block on CLL, Vc, Vp.
    # Table 2 reports %CV; the residual row shows the paper's convention
    # (sigma^2 = 0.0318 printed as %CV = 17.8 = 100 * sqrt(0.0318)), so
    # omega^2 = (%CV / 100)^2, the same convention as the predecessor
    # Rosario 2015 analysis by the same group.
    #   var(CLL) = 0.308^2 = 0.094864; var(Vc) = 0.202^2 = 0.040804;
    #   var(Vp)  = 0.702^2 = 0.492804
    #   cov(CLL, Vc) = 0.581  * 0.308 * 0.202 = 0.036148
    #   cov(CLL, Vp) = 0.0188 * 0.308 * 0.702 = 0.004065
    #   cov(Vc, Vp)  = 0.371  * 0.202 * 0.702 = 0.052609
    # IIV on Vmax, Q and Km fixed to 0 (Table 2) -- not encoded as etas.
    etalcl + etalvc + etalvp ~ c(
      0.094864,
      0.036148, 0.040804,
      0.004065, 0.052609, 0.492804
    ) # Table 2: IIV CLL 30.8 %CV, IIV Vc 20.2 %CV, IIV Vp 70.2 %CV; CORR CLL-Vc r = 0.581, CLL-Vp r = 0.0188, Vc-Vp r = 0.371

    # Inter-occasion variability on CLL -- Table 2: 20.3 %CV, omega^2 = 0.203^2
    # = 0.041209. One eta per occasion sharing the variance (NONMEM
    # $OMEGA BLOCK(1) SAME); occasions 2-10 are fixed to the occasion-1 value.
    etaiov_cl_1 ~ 0.041209 # Table 2: IOV CLL 20.3 %CV
    etaiov_cl_2 ~ fixed(0.041209) # same variance as occasion 1
    etaiov_cl_3 ~ fixed(0.041209) # same variance as occasion 1
    etaiov_cl_4 ~ fixed(0.041209) # same variance as occasion 1
    etaiov_cl_5 ~ fixed(0.041209) # same variance as occasion 1
    etaiov_cl_6 ~ fixed(0.041209) # same variance as occasion 1
    etaiov_cl_7 ~ fixed(0.041209) # same variance as occasion 1
    etaiov_cl_8 ~ fixed(0.041209) # same variance as occasion 1
    etaiov_cl_9 ~ fixed(0.041209) # same variance as occasion 1
    etaiov_cl_10 ~ fixed(0.041209) # same variance as occasion 1

    # Residual error -- Table 2: sigma^2_prop = 0.0318 (%CV = 17.8).
    propSd <- sqrt(0.0318); label("Proportional residual error on vedolizumab concentration (fraction)") # Table 2: Res prop sigma^2 = 0.0318 (%CV = 17.8)
  })

  model({
    # 1. Derived covariate terms.
    # AVA status from the time-varying titer: positive when titer >= 10
    # (Results 'Final PK Model'). AVA-negative records take the reference
    # titer 250 in the power term so it collapses to 1 and a 0-coded titer
    # cannot produce 0^e.
    ada_pos <- 0
    if (ADA_TITER >= 10) ada_pos <- 1
    ada_titer_use <- ada_pos * ADA_TITER + (1 - ada_pos) * 250

    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    oc5 <- (OCC == 5)
    oc6 <- (OCC == 6)
    oc7 <- (OCC == 7)
    oc8 <- (OCC == 8)
    oc9 <- (OCC == 9)
    oc10 <- (OCC == 10)
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 + oc3 * etaiov_cl_3 +
      oc4 * etaiov_cl_4 + oc5 * etaiov_cl_5 + oc6 * etaiov_cl_6 +
      oc7 * etaiov_cl_7 + oc8 * etaiov_cl_8 + oc9 * etaiov_cl_9 +
      oc10 * etaiov_cl_10

    # 2. Individual parameters. Separate typical CLL for AVA-negative and
    # AVA-positive records sharing the IIV, IOV and the WT / albumin / race /
    # diagnosis effects (Table 2 footnote a). Albumin reference 4 g/dL = 40 g/L.
    cl_typ <- (1 - ada_pos) * exp(lcl) +
      ada_pos * exp(lcl_adapos) * (ada_titer_use / 250)^e_ada_titer_cl
    cl <- cl_typ * exp(etalcl + iov_cl) *
      (WT / 70)^e_wt_cl *
      (ALB / 40)^e_alb_cl *
      e_race_asian_cl^RACE_ASIAN *
      e_ibd_cd_cl^IBD_CD
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vp
    q <- exp(lq) * (WT / 70)^e_wt_q
    vmax <- exp(lvmax) * (WT / 70)^e_wt_vmax
    km <- exp(lkm)

    # 3. Micro-constants.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. Two-compartment model, zero-order IV infusion into the central
    # compartment, parallel linear and Michaelis-Menten elimination
    # (Okamoto 2021 Figure 1; ADVAN13). Dose mg, volume L -> mg/L = ug/mL.
    Cc <- central / vc
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1 -
      vmax * Cc / (km + Cc)
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # 6. Observation and error.
    Cc ~ prop(propSd)
  })
}
