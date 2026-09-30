Schmidt_2020_camptothecin <- function() {
  description <- paste(
    "Coupled two-analyte population PK model of nanoparticle-bound",
    "(conjugated; output Cc_np) and free (released; output Cc) camptothecin",
    "after intravenous NLG207 (formerly CRLX101), a cyclodextrin-polymer",
    "nanoparticle-drug conjugate of camptothecin, in adults with advanced",
    "solid tumours (Schmidt 2020). Conjugated camptothecin follows a linear",
    "two-compartment model whose only elimination is conversion to free",
    "camptothecin in the central compartment, at a release clearance that",
    "is the sum of a steady-state arm and a fast arm decaying",
    "mono-exponentially with time after dose (half-life 0.307 h). Free",
    "camptothecin follows a second linear two-compartment model with",
    "first-order elimination. All clearances and volumes except the fast",
    "release arm and its decay are allometrically scaled on body weight",
    "with fixed exponents (0.75 and 1, reference 70 kg).",
    sep = " "
  )
  reference <- paste(
    "Schmidt KT, Huitema ADR, Dorlo TPC, Peer CJ, Cordes LM, Sciuto L,",
    "Wroblewski S, Pommier Y, Madan RA, Thomas A, Figg WD. Population",
    "pharmacokinetic analysis of nanoparticle-bound and free camptothecin",
    "after administration of NLG207 in adults with advanced solid tumors.",
    "Cancer Chemother Pharmacol. 2020;86(4):475-486.",
    "doi:10.1007/s00280-020-04134-9. PMCID: PMC7515962.",
    sep = " "
  )
  vignette <- "Schmidt_2020_camptothecin"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Actual body weight; allometric scaling with reference 70 kg and",
        "fixed exponents 1 (all four volumes) and 0.75 (Q1, CL_B, Q3, CL3)",
        "per Schmidt 2020 Table 2 parameter formulas. The fast release arm",
        "CL_F and its decay half-life are not weight-scaled (Results,",
        "'Covariate model'). Cohort median 70.4 kg (range 46.4-105 kg,",
        "Table 1)."
      ),
      source_name = "BW"
    )
  )

  compartmentData <- list(
    central_np = list(
      analyte = "camptothecin (nanoparticle-conjugated, NLG207)",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1_np = list(
      analyte = "camptothecin (nanoparticle-conjugated, NLG207)",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    central = list(analyte = "camptothecin (free)", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "camptothecin (free)", units = "mg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 27L,
    n_studies = 2L,
    n_observations = 477L,
    age_range = "47-76 years",
    age_median = "60 years",
    weight_range = "46.4-105 kg",
    weight_median = "70.4 kg",
    sex_female_pct = 55.6,
    race_ethnicity = c(Caucasian = 74.1, African = 14.8, Asian = 11.1),
    disease_state = paste(
      "Advanced solid tumours (NSCLC, small cell, pancreatic,",
      "cholangiocarcinoma, ovarian/fallopian tube, mCRPC, cervical,",
      "colorectal, mesothelioma, myxofibrosarcoma, thymic)"
    ),
    dose_range = paste(
      "NLG207 12 mg/m^2 IV over 1 or 2 h every 2 weeks at cycle 1",
      "(two cycle-6 patients dose-reduced to 50% and 75%)"
    ),
    renal_function = "eGFR 60-90 mL/min/1.73 m^2 in 48.1%; > 90 in 51.9%",
    co_medication = paste(
      "Olaparib (NCT02769962, started 48 h after the NLG207 infusion) or",
      "enzalutamide (NCT03531827); cycle-1 PK samples were collected without",
      "the combination partner"
    ),
    regions = "United States (National Cancer Institute)",
    notes = paste(
      "Schmidt 2020 Table 1. Two phase II NCI trials (NCT02769962 NLG207 +",
      "olaparib, n = 24; NCT03531827 NLG207 + enzalutamide, n = 3). 239",
      "plasma samples over 32 doses, each assayed for conjugated and free",
      "camptothecin; five patients also sampled at cycle 6."
    )
  )

  ini({
    # Conjugated (nanoparticle-bound) camptothecin -- Schmidt 2020 Table 2,
    # typical values for a 70 kg subject.
    lvc_np <- log(3.16); label("Central volume of conjugated CPT, V1 (L)")                          # Table 2 'V1,70kg' = 3.16 L
    lvp_np <- log(2.09); label("Peripheral volume of conjugated CPT, V2 (L)")                       # Table 2 'V2,70kg' = 2.09 L
    lq_np <- log(0.0482); label("Intercompartmental clearance of conjugated CPT, Q1 (L/h)")        # Table 2 'Q1,70kg' = 0.0482 L/h
    lcl_exp_inf_np <- log(0.0988); label("Steady-state (base) CPT release clearance, CL_B (L/h)")  # Table 2 'CLB,70kg' = 0.0988 L/h
    lcl_exp_component_np <- log(5.71); label("Initial fast CPT release clearance, CL_F (L/h)")     # Table 2 'CLF' = 5.71 L/h
    # Table 2 reports the decay as a half-life t1/2 = 0.307 h; the rate
    # constant is ln(2) / 0.307 = 2.2578 1/h (Fig. 1 CL1 equation).
    lcl_exp_kdes_np <- log(log(2) / 0.307); label("Decay rate constant of the fast release arm, ln(2)/t1/2 (1/h)") # Table 2 't1/2' = 0.307 h

    # Free (released) camptothecin -- Schmidt 2020 Table 2 (70 kg). Estimated
    # relative to the unknown fraction of conjugated CPT converted (Methods).
    lvc <- log(21.1); label("Central volume of free CPT, V3 (L)")                                  # Table 2 'V3,70kg' = 21.1 L
    lvp <- log(19.4); label("Peripheral volume of free CPT, V4 (L)")                               # Table 2 'V4,70kg' = 19.4 L
    lq <- log(25.6); label("Intercompartmental clearance of free CPT, Q3 (L/h)")                   # Table 2 'Q3,70kg' = 25.6 L/h
    lcl <- log(0.874); label("Clearance of free CPT from V3, CL3 (L/h)")                           # Table 2 'CL3,70kg' = 0.874 L/h

    # Allometric exponents fixed to standard values (Table 2 formulas; Results
    # 'Covariate model': estimated exponents were ~1.00 and reverted).
    e_wt_cl_q <- fixed(0.75); label("Allometric exponent on Q1, CL_B, Q3, CL3 (unitless)")        # Table 2 formulas '(BW/70kg)^0.75'
    e_wt_vc_vp <- fixed(1); label("Allometric exponent on V1, V2, V3, V4 (unitless)")               # Table 2 formulas '(BW/70kg)^1'

    # IIV: omega^2 = log(CV^2 + 1); covariance = r * sqrt(omega1^2 * omega2^2).
    # V1 CV 18.1% -> 0.0322358; CL_B CV 33.5% -> 0.106363; r = 0.918 -> 0.0537534
    etalvc_np + etalcl_exp_inf_np ~ c(0.0322358, 0.0537534, 0.106363)   # Table 2 'BSV V1' 18.1%, 'BSV V1,CLB (corr.)' 0.918, 'BSV CLB' 33.5%
    etalcl_exp_component_np ~ 0.330652                                   # Table 2 'BSV CLF' 62.6% CV
    # V3 CV 79.8% -> 0.492746; CL3 CV 42.2% -> 0.163889; r = 0.884 -> 0.251211
    etalvc + etalcl ~ c(0.492746, 0.251211, 0.163889)                    # Table 2 'BSV V3' 79.8%, 'BSV V3,CL3 (corr.)' 0.884, 'BSV CL3' 42.2%

    # Residual error: Eq. 2, Cobs = Cpred * (1 + eps_prop) + eps_add, per analyte.
    propSd_np <- 0.123; label("Proportional residual error, conjugated CPT (fraction)")            # Table 2 'BoundRE proportional' 12.3%
    addSd_np <- 5.07; label("Additive residual error, conjugated CPT (ng/mL)")                      # Table 2 'BoundRE additive' 5.07 ng/mL
    propSd <- 0.248; label("Proportional residual error, free CPT (fraction)")                      # Table 2 'FreeRE proportional' 24.8%
    addSd <- 0.396; label("Additive residual error, free CPT (ng/mL)")                              # Table 2 'FreeRE additive' 0.396 ng/mL
  })

  model({
    # Allometric size terms (reference 70 kg)
    wt_cl <- (WT / 70)^e_wt_cl_q
    wt_v <- (WT / 70)^e_wt_vc_vp

    # Conjugated CPT
    vc_np <- exp(lvc_np + etalvc_np) * wt_v
    vp_np <- exp(lvp_np) * wt_v
    q_np <- exp(lq_np) * wt_cl
    cl_exp_inf_np <- exp(lcl_exp_inf_np + etalcl_exp_inf_np) * wt_cl
    cl_exp_component_np <- exp(lcl_exp_component_np + etalcl_exp_component_np)
    cl_exp_kdes_np <- exp(lcl_exp_kdes_np)

    # Free CPT
    vc <- exp(lvc + etalvc) * wt_v
    vp <- exp(lvp) * wt_v
    q <- exp(lq) * wt_cl
    cl <- exp(lcl + etalcl) * wt_cl

    # Fig. 1: CL1 = CL_B + CL_F * exp(-(ln 2 / t1/2) * t), with t the time since
    # the start of the NLG207 infusion (Discussion: the release half-life is
    # 0.38 h "immediately at start of infusion" and 22 h by 4 h post-start).
    tdose <- tad()
    cl_release <- cl_exp_inf_np + cl_exp_component_np * exp(-cl_exp_kdes_np * tdose)

    # ODE system; release of CPT occurs only from the conjugated central
    # compartment into the free central compartment (Methods).
    d/dt(central_np) <- -cl_release / vc_np * central_np - q_np / vc_np * central_np + q_np / vp_np * peripheral1_np
    d/dt(peripheral1_np) <- q_np / vc_np * central_np - q_np / vp_np * peripheral1_np
    d/dt(central) <- cl_release / vc_np * central_np - cl / vc * central - q / vc * central + q / vp * peripheral1
    d/dt(peripheral1) <- q / vc * central - q / vp * peripheral1

    # Observations: dose in mg, volume in L -> mg/L; x 1000 -> ng/mL
    Cc_np <- 1000 * central_np / vc_np
    Cc <- 1000 * central / vc

    Cc_np ~ add(addSd_np) + prop(propSd_np)
    Cc ~ add(addSd) + prop(propSd)
  })
}
