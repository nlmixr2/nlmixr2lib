Kleideiter_2018_cebranopadol <- function() {
  description <- "Two-compartment population PK model for cebranopadol, a NOP / opioid receptor agonist, in healthy adults and adult chronic-pain patients (low back pain or osteoarthritis, diabetic polyneuropathy, post-bunionectomy), with two transit absorption compartments before central, first-order elimination, and covariate effects from sex, CYP2C9 phenotype, ALT, CrCl, age, body weight, formulation, and disease status (Kleideiter 2018)"
  reference   <- paste(
    "Kleideiter E, Piana C, Wang S, Nemeth R, Gautrois M. Clinical pharmacokinetic characteristics of cebranopadol, a novel first-in-class analgesic. Clin Pharmacokinet. 2018;57(1):31-50. doi:10.1007/s40262-017-0545-1.",
    "Erratum in: Kleideiter E, Piana C, Wang S, Nemeth R, Gautrois M. Correction to: clinical pharmacokinetic characteristics of cebranopadol, a novel first-in-class analgesic. Clin Pharmacokinet. 2018;57(11):1467-1469. doi:10.1007/s40262-018-0686-x.",
    sep = " "
  )
  vignette    <- "Kleideiter_2018_cebranopadol"
  units       <- list(time = "hr", dosing = "ug", concentration = "pg/mL")

  covariateData <- list(
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on Vc with reference 55 years (Kleideiter 2018 Table 14 / erratum reference covariate values). Estimated exponent -0.446 (Table 13).",
      source_name        = "AGE"
    ),
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on Vp with reference 82 kg (Kleideiter 2018 Table 14 / erratum reference covariate values). Estimated exponent 0.604 (Table 13). The paper tested but did not retain a body-weight effect on CL.",
      source_name        = "WT"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (female; the most-common category in the pooled analysis cohort and the Table 13 typical-value reference for CL).",
      notes              = "Multiplicative effect on CL applied as ratio^(1 - SEXF) so that female (SEXF = 1) is the reference (factor 1) and male (SEXF = 0) carries the +17.6% male-vs-female CL multiplier (87.4 / 74.3 = 1.176; Kleideiter 2018 Table 13).",
      source_name        = "SEX"
    ),
    CRCL = list(
      description        = "Cockcroft-Gault creatinine clearance",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL with reference 106.4 mL/min (Kleideiter 2018 Table 14 / erratum reference covariate values). Estimated exponent 0.349 (Table 13); higher CrCl is associated with higher cebranopadol CL.",
      source_name        = "CRCL"
    ),
    ALT = list(
      description        = "Baseline alanine aminotransferase activity",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL with reference 19 U/L (Kleideiter 2018 Table 14 / erratum reference covariate values). Estimated exponent -0.156 (Table 13); higher ALT is associated with lower cebranopadol CL.",
      source_name        = "ALT"
    ),
    FORM_CAPSULE = list(
      description        = "Liquid-filled capsule formulation indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (the typical-value reference is the film-coated tablet, with the sibling FORM_SOLUTION canonical also at 0 in the reference state).",
      notes              = "Multiplicative effects: on Ka (ratio 2.09 / 0.864 = 2.42), on klag (ratio 0.077 / 0.087 = 0.885), and on bioavailability F (factor 1.174 vs the tablet reference); Kleideiter 2018 Table 13.",
      source_name        = "FORM"
    ),
    FORM_SOLUTION = list(
      description        = "Oral solution formulation indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (the typical-value reference is the film-coated tablet, with the sibling FORM_CAPSULE canonical also at 0 in the reference state).",
      notes              = "Multiplicative effects: on Ka (ratio 2.43 / 0.864 = 2.81), on klag (ratio 0.077 / 0.087 = 0.885), and on bioavailability F (factor 1.045 vs the tablet reference); Kleideiter 2018 Table 13.",
      source_name        = "FORM"
    ),
    CYP2C9_EM = list(
      description        = "CYP2C9 extensive-metabolizer phenotype indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (in the Kleideiter 2018 pooled cohort the 0-level used here pools CYP2C9 phenotype unknown together with non-EM phenotypes, matching the model's Table 13 reference 'CYP2C9 status = unknown'; differs from the Jeong 2022 use of CYP2C9_EM where 0 = IM or PM only).",
      notes              = "Only 38.3% of the analysis cohort had a known CYP2C9 phenotype; the remaining 61.7% were assigned 'unknown' which is the model's reference category. Multiplicative effect on CL: 82.4 / 74.3 = 1.109 (about +11% CL vs unknown reference); Kleideiter 2018 Table 13.",
      source_name        = "CYP2C9"
    ),
    CYP2C9_PM_IM = list(
      description        = "CYP2C9 poor- or intermediate-metabolizer phenotype indicator (pooled)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (in the Kleideiter 2018 pooled cohort the 0-level pools CYP2C9 phenotype unknown together with EM, matching the model's Table 13 reference 'CYP2C9 status = unknown').",
      notes              = "Pooled poor + intermediate metabolizer indicator; the paper does not separately resolve PM from IM. Multiplicative effect on CL: 58.7 / 74.3 = 0.790 (about -21% CL vs unknown reference); Kleideiter 2018 Table 13.",
      source_name        = "CYP2C9"
    ),
    DIS_HEALTHY = list(
      description        = "Healthy-participant cohort indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (in this model 0 = the pooled nociceptive-pain reference category, low back pain or osteoarthritis, matching Kleideiter 2018 Table 13 / Table 14 erratum 'disease status = LBP and OA patients').",
      notes              = "Multiplicative effect on F: factor 0.837 (about -16% bioavailability vs the OA/LBP reference); Kleideiter 2018 Table 13.",
      source_name        = "DISEASE"
    ),
    DIS_DPN = list(
      description        = "Diabetic polyneuropathy patient indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (the pooled nociceptive-pain reference category, low back pain or osteoarthritis; matches the Kleideiter 2018 Table 13 / Table 14 erratum 'disease status = LBP and OA patients').",
      notes              = "Multiplicative effect on F: factor 1.132 (about +13% bioavailability vs the OA/LBP reference); Kleideiter 2018 Table 13 (erratum-corrected: the original Table 13 row labels for DPN and Bunionectomy were swapped per the Page 46 correction in the 2018 erratum).",
      source_name        = "DISEASE"
    ),
    DIS_BUN = list(
      description        = "Post-bunionectomy patient indicator (first-metatarsal bunionectomy cohort)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (the pooled nociceptive-pain reference category, low back pain or osteoarthritis; matches the Kleideiter 2018 Table 13 / Table 14 erratum 'disease status = LBP and OA patients').",
      notes              = "Multiplicative effect on F: factor 1.801 (about +80% bioavailability vs the OA/LBP reference); Kleideiter 2018 Table 13 (erratum-corrected: the original Table 13 row labels for Bunionectomy and DPN were swapped per the Page 46 correction in the 2018 erratum).",
      source_name        = "DISEASE"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 1293L,                            # 287 phase I + 1006 phase II; Kleideiter 2018 Results / Table 12 pooled across 14 trials
    n_studies      = 14L,                              # 8 phase I + 6 phase II; Kleideiter 2018 Tables 1-2
    age_range      = "18-79 years (median 33 in phase I, 58 in phase II)",  # Kleideiter 2018 Results 3.2.7 paragraph 1
    age_median     = "Phase I 33 y, Phase II 58 y",
    weight_range   = "approximately 45-197 kg",        # Kleideiter 2018 Table 12 min/max across all trials
    weight_median  = "approximately 82 kg",            # Kleideiter 2018 Table 14 erratum reference covariate
    sex_female_pct = 52.3,                             # 676 women / 1293 total; Kleideiter 2018 Discussion paragraph 4 (617 men, 676 women)
    race_ethnicity = "Not reported in detail; the analysis cohort is the union of European, North American, and other phase I / phase II trial enrollment.",
    disease_state  = "Pooled healthy adults (phase I) and adult chronic-pain patients across four indications: chronic low back pain (cLBP), osteoarthritis (OA) of the knee, painful diabetic polyneuropathy (DPN), and post-bunionectomy acute pain. The Table 13 typical-value reference for disease status is the nociceptive-pain group (OA + cLBP).",
    dose_range     = "Single oral doses 0.8-800 ug across multiple formulations (oral solution, film-coated tablet, liquid-filled capsule, encapsulated tablet); multiple oral doses 25-1600 ug once daily for 5 days to 14 weeks. The simulated typical therapeutic dose is 600 ug/day reached by a 100 -> 200 -> 400 -> 600 ug titration scheme (Kleideiter 2018 Table 14 footnote, erratum-corrected).",
    regions        = "Multinational phase I and II trials.",
    cyp2c9_status  = "Genotype-derived CYP2C9 phenotype available for 38.3% of subjects (extensive or poor / intermediate metabolizer); the remaining 61.7% are 'unknown' which is the model's typical-value reference for CYP2C9.",
    notes          = "Population characteristics summarised from Kleideiter 2018 Table 12 (per-trial demographics), Section 3.2.7 paragraph 1 (combined phase-I / phase-II counts and age ranges), and the Discussion paragraph 4 (sex split). Reference covariate values from the Table 14 footnote (erratum-corrected): sex = male, formulation = tablet, CYP2C9 status = unknown, disease status = LBP and OA patients, age (years) = 55, CrCl (mL/min) = 106.4, body weight (kg) = 82, ALT (units/L) = 19. Note that the Table 13 typical-value reference for sex is female (the most-common category), while the Table 14 simulation reference is the male typical patient; the model file uses the Table 13 (parameter-estimation) reference, with male encoded as a (1 - SEXF) effect."
  )

  ini({
    # Structural PK parameters at the Table 13 reference covariate set:
    #   SEXF = 1 (female; the Table 13 most-common-category reference), AGE = 55 y,
    #   WT = 82 kg, CRCL = 106.4 mL/min, ALT = 19 U/L,
    #   FORM_CAPSULE = FORM_SOLUTION = 0 (tablet), CYP2C9_EM = CYP2C9_PM_IM = 0 (unknown),
    #   DIS_HEALTHY = DIS_DPN = DIS_BUN = 0 (LBP / OA reference).
    lcl    <- log(74.3);     label("Apparent oral clearance at the reference covariate set (CL/F, L/h)")                       # Kleideiter 2018 Table 13 'Clearance Reference value'
    lvc    <- log(225);      label("Apparent central volume of distribution at the reference covariate set (Vc/F, L)")         # Kleideiter 2018 Table 13 'Volume central compartment Reference value'
    lvp    <- log(6750);     label("Apparent peripheral volume of distribution at the reference covariate set (Vp/F, L)")      # Kleideiter 2018 Table 13 'Volume peripheral compartment Reference value'
    lq     <- log(84.2);     label("Apparent inter-compartmental clearance at the reference covariate set (Q/F, L/h)")         # Kleideiter 2018 Table 13 'Intercompartmental clearance'
    lka    <- log(0.864);    label("Absorption rate constant at the reference covariate set (Ka, 1/h)")                        # Kleideiter 2018 Table 13 'Absorption rate constant Reference value'
    lklag  <- log(0.087);    label("Transition-compartment rate constant at the reference covariate set (klag, 1/h)")          # Kleideiter 2018 Table 13 'klag Reference value'

    # Covariate effects on CL (Kleideiter 2018 Table 13)
    e_male_cl   <- 1.176;    label("Male-vs-female CL ratio (87.4 / 74.3); applied as ratio^(1 - SEXF)")                       # Kleideiter 2018 Table 13: Male CL 87.4 vs Female CL 74.3 = 1.176
    e_em_cl     <- 1.109;    label("CYP2C9 EM-vs-unknown CL ratio (82.4 / 74.3); applied as ratio^CYP2C9_EM")                  # Kleideiter 2018 Table 13: EM CL 82.4 vs Reference 74.3 = 1.109
    e_pim_cl    <- 0.790;    label("CYP2C9 PIM-vs-unknown CL ratio (58.7 / 74.3); applied as ratio^CYP2C9_PM_IM")                # Kleideiter 2018 Table 13: PIM CL 58.7 vs Reference 74.3 = 0.790
    e_alt_cl    <- -0.156;   label("Power exponent of ALT on CL (unitless; reference 19 U/L)")                                 # Kleideiter 2018 Table 13 'Effect of ALT (exponential)'
    e_crcl_cl   <-  0.349;   label("Power exponent of CrCl on CL (unitless; reference 106.4 mL/min)")                          # Kleideiter 2018 Table 13 'Effect of CrCl (exponential)'

    # Covariate effects on Vc (Kleideiter 2018 Table 13)
    e_age_vc    <- -0.446;   label("Power exponent of AGE on Vc (unitless; reference 55 years)")                               # Kleideiter 2018 Table 13 'Effect of age (exponential)'

    # Covariate effects on Vp (Kleideiter 2018 Table 13)
    e_wt_vp     <-  0.604;   label("Power exponent of WT on Vp (unitless; reference 82 kg)")                                   # Kleideiter 2018 Table 13 'Effect of body weight (exponential)'

    # Covariate effects on Ka (Kleideiter 2018 Table 13)
    e_solution_ka <- 2.813;  label("Oral-solution-vs-tablet Ka ratio (2.43 / 0.864); applied as ratio^FORM_SOLUTION")          # Kleideiter 2018 Table 13: Oral solution 2.43 vs Reference 0.864 = 2.813
    e_capsule_ka  <- 2.419;  label("Capsule-vs-tablet Ka ratio (2.09 / 0.864); applied as ratio^FORM_CAPSULE")                 # Kleideiter 2018 Table 13: Capsules 2.09 vs Reference 0.864 = 2.419

    # Covariate effects on klag (Kleideiter 2018 Table 13)
    e_solution_klag <- 0.885; label("Oral-solution-vs-tablet klag ratio (0.077 / 0.087); applied as ratio^FORM_SOLUTION")      # Kleideiter 2018 Table 13: Oral solution 0.077 vs Reference 0.087 = 0.885
    e_capsule_klag  <- 0.885; label("Capsule-vs-tablet klag ratio (0.077 / 0.087); applied as ratio^FORM_CAPSULE")             # Kleideiter 2018 Table 13: Capsules 0.077 vs Reference 0.087 = 0.885

    # Covariate effects on bioavailability F (Kleideiter 2018 Table 13)
    # The structural F is fixed at 1 at the reference covariate set; covariate-driven
    # multiplicative shifts express the formulation and disease-state effects.
    e_solution_f <- 1.045;   label("Oral-solution-vs-tablet F ratio; applied as ratio^FORM_SOLUTION")                          # Kleideiter 2018 Table 13 'Bioavailability Oral solution'
    e_capsule_f  <- 1.174;   label("Capsule-vs-tablet F ratio; applied as ratio^FORM_CAPSULE")                                 # Kleideiter 2018 Table 13 'Bioavailability Capsules'
    e_healthy_f  <- 0.837;   label("Healthy-vs-LBP/OA F ratio; applied as ratio^DIS_HEALTHY")                                  # Kleideiter 2018 Table 13 'Bioavailability Healthy volunteers'
    # Erratum Page 46 correction: the original Table 13 row labels for
    # 'Bunionectomy patients' (F = 1.132) and 'DPN patients' (F = 1.801) were
    # swapped; the corrected mapping places F = 1.801 with bunionectomy and
    # F = 1.132 with DPN. The erratum's Table 14 (Cmax,ss +80% / AUCss +80% for
    # bunionectomy; +13% / +13% for DPN) confirms this assignment. Per the
    # Phase 1 step 8 errata-handling convention, the erratum-corrected mapping
    # is the source of record here.
    e_bun_f      <- 1.801;   label("Bunionectomy-vs-LBP/OA F ratio; applied as ratio^DIS_BUN")                                 # Kleideiter 2018 Table 13 (erratum-corrected; original-Table-13-row-labeled 'DPN patients')
    e_dpn_f      <- 1.132;   label("DPN-vs-LBP/OA F ratio; applied as ratio^DIS_DPN")                                          # Kleideiter 2018 Table 13 (erratum-corrected; original-Table-13-row-labeled 'Bunionectomy patients')

    # Inter-individual variability. Kleideiter 2018 Table 13 reports log-normal omega^2
    # values directly (column 'Interindividual variability (RSE%)'); no off-diagonal
    # covariances are reported, so a diagonal omega matrix is used.
    etalcl    ~ 0.412     # Kleideiter 2018 Table 13: CL omega^2 = 0.412 (RSE 10.1%)
    etalvc    ~ 0.559     # Kleideiter 2018 Table 13: Vc omega^2 = 0.559 (RSE 20.6%)
    etalka    ~ 0.519     # Kleideiter 2018 Table 13: Ka omega^2 = 0.519 (RSE 11.2%)
    etalklag  ~ 0.0626    # Kleideiter 2018 Table 13: klag omega^2 = 0.0626 (RSE 18.8%)

    # Residual error. The paper (Kleideiter 2018 Section 2.4.1, Eq. 2) reports the
    # error model as additive on log-transformed cebranopadol concentrations (i.e.,
    # proportional on the linear scale) but does NOT report the residual-error
    # magnitude in Table 13 or the Results narrative. Approximated here from typical
    # popPK class values for oral small-molecule analgesics with mixed intensive +
    # sparse sampling; see vignette Errata / Assumptions and deviations.
    # paper does not report a residual-error magnitude; class-typical approximation
    propSd <- 0.35;       label("Proportional residual error (fraction; paper does not report, class-typical approximation)")
  })

  model({
    # Individual PK parameters at the per-subject covariate values, evaluated against
    # the Table 13 reference set (Section 2.4.1 Eq. 3 power-form covariate model for
    # continuous covariates; Eq. 4 linear-deviation form for categorical effects,
    # reformulated here as multiplicative `ratio^indicator` to keep the typical-value
    # ratios from Table 13 reproducible at a glance).

    cl <- exp(lcl + etalcl) *
      e_male_cl^(1 - SEXF) *
      e_em_cl^CYP2C9_EM *
      e_pim_cl^CYP2C9_PM_IM *
      (ALT  / 19   )^e_alt_cl *
      (CRCL / 106.4)^e_crcl_cl
    vc <- exp(lvc + etalvc) *
      (AGE / 55)^e_age_vc
    vp <- exp(lvp) *
      (WT / 82)^e_wt_vp
    q  <- exp(lq)
    ka <- exp(lka + etalka) *
      e_solution_ka^FORM_SOLUTION *
      e_capsule_ka^FORM_CAPSULE
    klag <- exp(lklag + etalklag) *
      e_solution_klag^FORM_SOLUTION *
      e_capsule_klag^FORM_CAPSULE
    fdepot <- e_solution_f^FORM_SOLUTION *
      e_capsule_f^FORM_CAPSULE *
      e_healthy_f^DIS_HEALTHY *
      e_bun_f^DIS_BUN *
      e_dpn_f^DIS_DPN

    # Micro-constants for the explicit ODE system.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ODE system. Structural interpretation of the absorption chain (Kleideiter 2018
    # Section 3.2.7: "two lagged transition compartments"): the dose enters the depot,
    # the formulation-dependent first-order rate Ka governs the depot -> transit1 step
    # (dissolution-limited absorption), and the transit-compartment rate klag governs
    # both subsequent first-order steps transit1 -> transit2 -> central. See vignette
    # Assumptions and deviations for the rationale and for the alternative reading
    # (klag at depot + transit1, Ka at transit2 -> central) that the paper text could
    # also support.
    d/dt(depot)       <- -ka * depot
    d/dt(transit1)      <-  ka  * depot  - klag * transit1
    d/dt(transit2)      <-  klag * transit1 - klag * transit2
    d/dt(central)     <-  klag * transit2 - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    f(depot) <- fdepot

    # Concentration: dose in ug, volume in L gives ug/L. Multiply by 1000 to express
    # Cc in pg/mL to match the units used throughout Kleideiter 2018 Tables 4-14.
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
