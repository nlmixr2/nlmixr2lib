Kleideiter_2017_cebranopadol <- function() {
  description <- "Two-compartment population PK model for oral cebranopadol with two lagged transition compartments in healthy subjects and chronic-pain patients (Kleideiter 2017; with 2018 correction)"
  reference <- paste(
    "Kleideiter E, Piana C, Wang S, Nemeth R, Gautrois M.",
    "Clinical Pharmacokinetic Characteristics of Cebranopadol, a Novel First-in-Class Analgesic.",
    "Clin Pharmacokinet. 2018;57(1):31-50. doi:10.1007/s40262-017-0545-1.",
    "Correction: Clin Pharmacokinet. 2018;57(11):1471-1472. doi:10.1007/s40262-018-0686-x.",
    "This model already incorporates the 2018 erratum: the Table 13 bioavailability rows for",
    "bunionectomy and DPN patients were swapped in the correction, so the packaged model uses",
    "the corrected assignment F_bunionectomy = 1.801 and F_DPN = 1.132."
  )
  vignette <- "Kleideiter_2017_cebranopadol"
  units <- list(time = "hour", dosing = "ug", concentration = "pg/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on Vp/F with reference weight 82 kg (Table 14 reference covariate values).",
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on Vc/F with reference age 55 years (Table 14 reference covariate values).",
      source_name        = "AGE"
    ),
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (female; the paper's reference category was female and the canonical reference is male, so the model carries a derived `male` indicator and applies the male effect when SEXF = 0)",
      notes              = "Kleideiter 2017 uses female sex as the typical-value reference (CL_ref = 74.3 L/h) and reports a male CL of 87.4 L/h. The model derives a male indicator as (1 - SEXF) and applies `e_male_lcl = log(87.4 / 74.3) = 0.162` to it.",
      source_name        = "SEX"
    ),
    CRCL = list(
      description        = "Cockcroft-Gault creatinine clearance (NOT BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Methods Section 2.4.1 states CrCl was derived from observed data via Cockcroft-Gault (raw mL/min). Power scaling on CL with reference 106.4 mL/min (Table 14 reference covariate value). Distinct from the canonical CRCL unit (mL/min/1.73 m^2) - follow the `Delattre_2010_amikacin` precedent of using `CRCL` with a per-model unit note.",
      source_name        = "CRCL"
    ),
    ALT = list(
      description        = "Alanine aminotransferase",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on CL with reference 19 U/L (Table 14 reference covariate value).",
      source_name        = "ALT"
    ),
    CYP2C9_EM = list(
      description        = "CYP2C9 extensive-metabolizer phenotype indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (intermediate / poor metabolizer OR unknown phenotype, when paired with CYP2C9_PM_IM)",
      notes              = "Kleideiter 2017 has three CYP2C9 strata: unknown (most common, reference), extensive metabolizer (EM), and pooled poor/intermediate metabolizer (PM/IM). CYP2C9_EM = 1 only for confirmed EM; CYP2C9_PM_IM = 1 only for confirmed PM/IM; both = 0 indicates unknown phenotype (which the paper's covariate model treats as the typical-value reference, CL = 74.3 L/h).",
      source_name        = "CYP2C9"
    ),
    CYP2C9_PM_IM = list(
      description        = "CYP2C9 poor-or-intermediate-metabolizer pooled phenotype indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (extensive metabolizer OR unknown phenotype)",
      notes              = "Pairs with CYP2C9_EM to encode the three-level Kleideiter 2017 CYP2C9 stratification (unknown reference / EM / PM-IM). PM and IM are pooled in Kleideiter; downstream papers that distinguish PM from IM should register a separate CYP2C9_PM canonical and split this group.",
      source_name        = "CYP2C9"
    ),
    FORM_TABLET = list(
      description        = "Film-coated tablet formulation indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (tablet is the typical-value reference for ka, klag, and bioavailability)",
      notes              = "Kleideiter 2017 has three formulation strata: tablet (reference), oral solution, and liquid-filled capsules. FORM_TABLET = 1 for tablet, FORM_CAPSULE = 1 for capsule, both = 0 indicates oral solution. Tablet is the reference category; absorption-rate, klag, and bioavailability covariate effects are zero at FORM_TABLET = 1, FORM_CAPSULE = 0.",
      source_name        = "FORMULATION"
    ),
    FORM_CAPSULE = list(
      description        = "Liquid-filled capsule formulation indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (the paired reference is tablet via FORM_TABLET = 1; the third level - oral solution - is encoded as FORM_TABLET = 0 AND FORM_CAPSULE = 0)",
      notes              = "Pairs with FORM_TABLET to encode the three-level Kleideiter 2017 formulation stratification.",
      source_name        = "FORMULATION"
    ),
    DIS_HEALTHY = list(
      description        = "Healthy participant indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (chronic-low-back-pain / osteoarthritis patient; the most common 'nociceptive pain' category and the bioavailability reference in Table 13)",
      notes              = "Pairs with DIS_DPN and DIS_BUNIONECTOMY to encode the four-level disease-status stratification (LBP/OA reference, healthy, DPN, bunionectomy).",
      source_name        = "DIS"
    ),
    DIS_DPN = list(
      description        = "Diabetic polyneuropathy indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (LBP/OA reference; complement includes healthy, bunionectomy, LBP, OA)",
      notes              = "Pairs with DIS_HEALTHY and DIS_BUNIONECTOMY to encode the Kleideiter 2017 disease-status stratification.",
      source_name        = "DIS"
    ),
    DIS_BUNIONECTOMY = list(
      description        = "Post-bunionectomy acute-pain indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (LBP/OA reference; complement includes healthy, DPN, LBP, OA)",
      notes              = "Pairs with DIS_HEALTHY and DIS_DPN to encode the Kleideiter 2017 disease-status stratification. The 2018 correction swaps the Table 13 row labels for bunionectomy and DPN bioavailability: bunionectomy patients have F = 1.801 (this entry) and DPN patients have F = 1.132 (corrected; the uncorrected paper had these values swapped).",
      source_name        = "DIS"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 1293L,
    n_studies      = 14L,
    age_range      = "18-79 years (Phase I 18-64, Phase II 18-79)",
    age_median     = "33 years Phase I; 58 years Phase II (Section 3.2.7)",
    weight_range   = "45.4-197 kg (across pooled Trials 1-14, Table 12)",
    weight_median  = "73.9-99 kg per trial (Table 12); covariate-model reference 82 kg",
    sex_female_pct = 52.3,
    race_ethnicity = "Not reported in this paper.",
    disease_state  = "Healthy adults (Phase I, Trials 1-8) and patients with chronic pain conditions: bunionectomy acute pain (Trial 9), diabetic polyneuropathy (Trials 10, 12, 14), osteoarthritis (Trial 11), chronic low back pain (Trials 4, 13).",
    dose_range     = "0.8 ug single dose to 1600 ug/day multiple dose, oral immediate-release (tablet, liquid-filled capsule, oral solution).",
    regions        = "Not specified (sponsor Gruenenthal GmbH and Forest Research Institute; EudraCT and ClinicalTrials.gov registrations).",
    notes          = "Pooled population PK analysis across 8 Phase I and 6 Phase II trials. Table 12 lists per-trial demographics; reference covariate values (Table 14): age = 55 y, CrCl = 106.4 mL/min, weight = 82 kg, ALT = 19 U/L; most common categorical values: female sex, tablet formulation, nociceptive-pain disease status (LBP/OA), unknown CYP2C9 phenotype. Sex composition (~52% female) is approximate from Table 12 per-trial counts."
  )

  ini({
    # Structural parameters (typical-value tablet reference; Table 13)
    lka       <- log(0.864) ; label("Absorption rate constant ka_ref - depot -> transit chain (1/h, tablet reference)")  # Table 13 row "Absorption rate constant - Reference value"
    lklag     <- log(0.087) ; label("Transition compartment rate constant klag_ref - transit2 -> central (1/h, tablet reference)")  # Table 13 row "k_lag - Reference value"
    lcl       <- log(74.3)  ; label("Apparent clearance CL/F at the typical female / unknown-CYP2C9 / LBP-OA reference subject (L/h)")  # Table 13 row "Clearance - Reference value"
    lvc       <- log(225)   ; label("Apparent central volume of distribution Vc/F at reference age 55 y (L)")  # Table 13 row "Volume central compartment - Reference value"
    lvp       <- log(6750)  ; label("Apparent peripheral volume of distribution Vp/F at reference body weight 82 kg (L)")  # Table 13 row "Volume peripheral compartment - Reference value"
    lq        <- log(84.2)  ; label("Apparent inter-compartmental clearance Q/F (L/h)")  # Table 13 row "Intercompartmental clearance"

    # Covariate effects on log(CL) - categorical (additive on log scale)
    e_male_lcl       <- 0.1620 ; label("Male sex effect on log(CL/F) - additive log shift = log(87.4 / 74.3)")     # Table 13 "Males 87.4 L/h" vs reference 74.3
    e_cyp2c9em_lcl   <- 0.1037 ; label("CYP2C9 extensive-metabolizer effect on log(CL/F) - additive log shift = log(82.4 / 74.3)")  # Table 13 "CYP2C9 extensive metabolizers 82.4 L/h"
    e_cyp2c9pmim_lcl <- -0.2353; label("CYP2C9 poor/intermediate-metabolizer effect on log(CL/F) - additive log shift = log(58.7 / 74.3)")  # Table 13 "CYP2C9 poor and intermediate metabolizers 58.7 L/h"

    # Covariate exponents on CL (continuous, power form on (cov / reference))
    e_alt_lcl  <- -0.156 ; label("ALT exponent on CL/F (unitless; reference ALT = 19 U/L)")    # Table 13 "Effect of ALT (exponential)"
    e_crcl_lcl <-  0.349 ; label("CrCl exponent on CL/F (unitless; reference CrCl = 106.4 mL/min)") # Table 13 "Effect of CrCl (exponential)"

    # Covariate exponent on Vc (continuous, power form on (cov / reference))
    e_age_lvc  <- -0.446 ; label("Age exponent on Vc/F (unitless; reference age = 55 y)")     # Table 13 "Effect of age (exponential)"

    # Covariate exponent on Vp (continuous, power form on (cov / reference))
    e_wt_lvp   <-  0.604 ; label("Body weight exponent on Vp/F (unitless; reference weight = 82 kg)") # Table 13 "Effect of body weight (exponential)"

    # Formulation effects on absorption parameters (categorical; additive on log scale)
    # Tablet is the reference (effect = 0 when FORM_TABLET = 1, FORM_CAPSULE = 0).
    e_form_sol_lka    <- 1.0341 ; label("Oral-solution effect on log(ka) - additive log shift = log(2.43 / 0.864)")   # Table 13 "Oral solution 2.43 h^-1" (absorption rate)
    e_form_cap_lka    <- 0.8829 ; label("Capsule effect on log(ka) - additive log shift = log(2.09 / 0.864)")          # Table 13 "Capsules 2.09 h^-1"
    e_form_sol_lklag  <- -0.1221; label("Oral-solution effect on log(klag) - additive log shift = log(0.077 / 0.087)") # Table 13 "Oral solution 0.077 h^-1" (klag)
    e_form_cap_lklag  <- -0.1221; label("Capsule effect on log(klag) - additive log shift = log(0.077 / 0.087)")       # Table 13 "Capsules 0.077 h^-1"

    # Bioavailability multipliers (categorical; multiplicative on F).
    # Reference is tablet + LBP/OA patient (F = 1.0).
    e_form_sol_f      <- 1.045 ; label("Oral-solution bioavailability multiplier (relative to tablet)")    # Table 13 "Oral solution 1.045"
    e_form_cap_f      <- 1.174 ; label("Capsule bioavailability multiplier (relative to tablet)")          # Table 13 "Capsules 1.174"
    e_dis_healthy_f   <- 0.837 ; label("Healthy-volunteer bioavailability multiplier (relative to LBP/OA)") # Table 13 "Healthy volunteers 0.837"
    e_dis_dpn_f       <- 1.132 ; label("DPN-patient bioavailability multiplier (relative to LBP/OA, corrected)")          # Table 13 + 2018 Correction (Table 13 rows 27-28 swap): DPN 1.132
    e_dis_buni_f      <- 1.801 ; label("Bunionectomy-patient bioavailability multiplier (relative to LBP/OA, corrected)") # Table 13 + 2018 Correction (Table 13 rows 27-28 swap): bunionectomy 1.801

    # Inter-individual variability (variances on log-scale parameters; Table 13 IIV column)
    etalcl   ~ 0.412   # Table 13 IIV for Clearance (variance on log CL, RSE 10.1%)
    etalvc   ~ 0.559   # Table 13 IIV for Volume central compartment (RSE 20.6%)
    etalka   ~ 0.519   # Table 13 IIV for Absorption rate constant (RSE 11.2%)
    etalklag ~ 0.0626  # Table 13 IIV for k_lag (RSE 18.8%)

    # Residual error placeholder.
    # NOT reported in Kleideiter 2017: the paper states (Eq. 2) that the residual error
    # is additive on log-transformed concentration with variance sigma^2, but the
    # numerical sigma value is not given in Table 13 or anywhere in the main text or
    # 2018 correction. Set to 0.20 (~20% CV in linear space) as a documented placeholder;
    # see vignette Errata and tracking/operator_followups.md for the gap.
    propSd   <- 0.20 ; label("Proportional residual error (placeholder; sigma not reported in Kleideiter 2017 - see vignette Errata)")
  })

  model({
    # 1. Derived covariate terms

    # Sex: paper reference is female (SEXF = 1 -> e_male_lcl * 0 = 0).
    # The canonical SEXF reference is 0 = male, so the male-effect indicator is (1 - SEXF).
    male_indicator <- 1 - SEXF

    # Disease-status reference is LBP/OA (all three indicators = 0).
    dis_lbp_oa <- 1 - DIS_HEALTHY - DIS_DPN - DIS_BUNIONECTOMY

    # Formulation reference is tablet (FORM_TABLET = 1, FORM_CAPSULE = 0).
    # Oral solution is encoded as FORM_TABLET = 0, FORM_CAPSULE = 0.
    is_solution <- (1 - FORM_TABLET) * (1 - FORM_CAPSULE)

    # Bioavailability: multiplicative formulation x disease (Table 13 footnote
    # "For bioavailability the reference value was set to 1").
    f_form    <- FORM_TABLET * 1 + is_solution * e_form_sol_f + FORM_CAPSULE * e_form_cap_f
    f_disease <- dis_lbp_oa * 1 + DIS_HEALTHY * e_dis_healthy_f + DIS_DPN * e_dis_dpn_f + DIS_BUNIONECTOMY * e_dis_buni_f
    f_total   <- f_form * f_disease

    # 2. Individual PK parameters
    ka   <- exp(lka + e_form_sol_lka * is_solution + e_form_cap_lka * FORM_CAPSULE + etalka)
    klag <- exp(lklag + e_form_sol_lklag * is_solution + e_form_cap_lklag * FORM_CAPSULE + etalklag)
    cl   <- exp(lcl + e_male_lcl * male_indicator +
                 e_cyp2c9em_lcl * CYP2C9_EM + e_cyp2c9pmim_lcl * CYP2C9_PM_IM +
                 etalcl) *
            (ALT / 19)^e_alt_lcl *
            (CRCL / 106.4)^e_crcl_lcl
    vc   <- exp(lvc + etalvc) * (AGE / 55)^e_age_lvc
    vp   <- exp(lvp)          * (WT / 82)^e_wt_lvp
    q    <- exp(lq)

    # 3. Micro-constants
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # 4. ODE system - two lagged transition compartments + 2-cmt disposition.
    # Section 3.2.7 prose: "A two-compartment disposition model with two lagged transition
    # compartments and a first-order elimination process best describes cebranopadol data."
    # The paper does not write out the ODEs; the implemented arrangement places ka on the
    # depot -> transit1 -> transit2 chain (Erlang-style transit) and klag on the final
    # transit2 -> central step. See vignette "Assumptions and deviations" for the
    # rationale and an enumeration of alternative interpretations.
    d/dt(depot)       <- -ka * depot
    d/dt(transit1)    <-  ka * depot   - ka   * transit1
    d/dt(transit2)    <-  ka * transit1 - klag * transit2
    d/dt(central)     <-  klag * transit2 - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # 5. Bioavailability (paper bioavailability multipliers apply to the dose via depot F)
    f(depot) <- f_total

    # 6. Observation and error.
    # Dose is in ug, volumes are in L: central/vc gives ug/L = ng/mL.
    # The paper reports concentrations in pg/mL, so multiply by 1000 to match.
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
