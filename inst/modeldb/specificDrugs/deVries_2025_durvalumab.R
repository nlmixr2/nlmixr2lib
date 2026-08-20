deVries_2025_durvalumab <- function() {
  description <- "Two compartment PK model of durvalumab (anti-PD-L1) with parallel linear and Michaelis-Menten elimination in non-small cell lung cancer, used for TDM-based tailored dosing simulations (de Vries 2025)"
  reference <- paste(
    "de Vries F, Franssen EJF, Smit AAJ, Moes DJAR, van der Wekken AJ,",
    "Oude Munnink T, Hendrikx JJMA, Dumoulin DW, Koolen SLW, Kievit W,",
    "van den Heuvel MM, ter Heine R. TDM-Based Tailored Dosing of",
    "Durvalumab in Lung Cancer Patients: A Comprehensive Population",
    "Pharmacokinetic-Pharmacoeconomic Evaluation. Clin Pharmacokinet.",
    "2025;64:1507-1515. doi:10.1007/s40262-025-01555-8.",
    "Every parameter value below is transcribed verbatim from the NONMEM",
    "control stream printed in Supplementary Model Code (ESM 3) of that",
    "paper. Nothing was estimated by de Vries 2025; the control stream is a",
    "$SIMULATION ONLYSIM problem whose structural parameters are hardcoded",
    "literals in $PK, inherited from the semi-mechanistic durvalumab model",
    "of Baverel PG, Dubois VFS, Jin CY, Zheng Y, Song X, Jin X, et al.",
    "Population pharmacokinetics of durvalumab in cancer patients and",
    "association with longitudinal biomarkers of disease status. Clin",
    "Pharmacol Ther. 2018;103(4):631-642. doi:10.1002/cpt.982.",
    sep = " "
  )
  vignette <- "deVries_2025_durvalumab"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power allometric effects on CL (exponent 0.389) and Vc (exponent 0.406), both normalized to the Baverel 2018 reference weight of 69.8 kg. Baseline (time-fixed) in the de Vries 2025 simulation. The simulated cohort median is 64.0 kg (ESM 4 Table S1), which is NOT the allometric reference; do not conflate the two.",
      source_name        = "WT"
    ),
    SEXF = list(
      description        = "Sex (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Multiplicative effects for females on CL (0.857), Vc (0.835) and Vp (0.795). The ESM 3 control stream derives the indicator as FLASEX = 1 when the source SEX column equals 0, i.e. the source encodes SEX = 0 as female; SEXF = 1 - SEX under that encoding. The direction (female lowers CL and both volumes after weight allometry) is corroborated by the sibling durvalumab model Ogasawara_2020_durvalumab.R (female CL 0.791, Vc 0.790). See the vignette Assumptions section.",
      source_name        = "SEX"
    ),
    ALB = list(
      description        = "Serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Centred linear effect on CL: (1 - 0.035 * (ALB - 38)). Reported and used in SI g/L, so no g/dL conversion is applied. Centring constant 38 g/L is the Baverel 2018 reference; the simulated cohort median is 39.27 g/L (ESM 4 Table S1).",
      source_name        = "ALB"
    ),
    CRCL = list(
      description        = "Creatinine clearance (raw Cockcroft-Gault style, NOT BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Centred linear effect on CL: (1 + 0.00149 * (CRCL - 85.65)). ESM 4 Table S1 reports 'Creatinine clearance (mL/min)' with no BSA normalization, so the per-model unit is mL/min rather than the register's default mL/min/1.73 m^2 (same precedent as Delattre_2010_amikacin.R). Centring constant 85.65 mL/min is the Baverel 2018 reference; the simulated cohort median is 86.43 mL/min.",
      source_name        = "CRCL"
    ),
    ECOG_GE1 = list(
      description        = "Baseline ECOG performance status >= 1 (1 = ECOG 1 or worse, 0 = ECOG 0)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ECOG performance status 0)",
      notes              = "Power form on CL: 0.937^ECOG_GE1, i.e. a 6.3 percent lower CL for ECOG >= 1. The ESM 3 control stream generates the flag with a $MIX block (P(1) = 0.75 for ECOG 0, P(2) = 0.25 for ECOG 1); that mixture is a simulation device for drawing the covariate, not a mixture PK model, and the cohort proportions it produces are reproduced exactly in ESM 4 Table S1 (750 / 250). Only ECOG 0 and 1 occur, so the ordinal WHO_PS column collapses to this binary indicator.",
      source_name        = "FLAECOG"
    ),
    TUMSZ = list(
      description        = "Baseline tumor size (sum of target-lesion diameters)",
      units              = "mm",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Centred linear effect on CL: (1 + 0.00178 * (TUMSZ - 74.8)). Centring constant 74.8 mm is the Baverel 2018 reference; the simulated cohort median is 50.09 mm (ESM 4 Table S1), so the typical simulated patient sits below the centring value and the term reduces CL.",
      source_name        = "TUMORSIZE"
    ),
    SPDL1 = list(
      description        = "Baseline soluble PD-L1",
      units              = "pg/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Centred linear effect on the Michaelis-Menten Vmax only (not on linear CL): (1 + 0.00336 * (SPDL1 - 124.8)). Centring constant 124.8 pg/mL is the Baverel 2018 reference; the simulated cohort median is 138.34 pg/mL (ESM 4 Table S1).",
      source_name        = "SPDL1"
    )
  )

  covariatesDataExcluded <- list(
    ADA_POS = list(
      description = "Anti-drug-antibody positivity",
      units       = "(binary)",
      type        = "binary",
      notes       = "The ESM 3 control stream sets ADA = 0 ('NO ADA') for every simulated subject and then never references ADA in any CL / V / Vmax expression, and ESM 4 Table S1 reports anti-drug antibodies as 0 for the whole cohort. The covariate is therefore documented for provenance only and is deliberately absent from model()."
    )
  )

  compartmentData <- list(
    central     = list(analyte = "durvalumab", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "durvalumab", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    cohort_type    = "virtual (Monte-Carlo simulated, not an observed patient cohort)",
    n_subjects     = 1000,
    n_studies      = 0,
    age_range      = "50-80 years",
    age_median     = "63.14 years",
    weight_range   = "40.6-89.4 kg",
    weight_median  = "64.0 kg",
    sex_female_pct = 50.5,
    disease_state  = "non-small cell lung cancer (ECOG performance status 0 in 75 percent, 1 in 25 percent)",
    dose_range     = "10 mg/kg IV Q2W and 1500 mg IV Q4W (approved regimens), plus TDM-tailored doses (1000-1740 mg) and TDM-tailored intervals (4-8 weeks); all administered as 1 h infusions",
    regions        = "the Netherlands (simulated cohort intended to represent Dutch NSCLC practice)",
    notes          = "Baseline characteristics are the simulated virtual population of ESM 4 Table S1 (N = 1000), not an enrolled cohort: albumin median 39.27 g/L [27.2-55.3], creatinine clearance median 86.43 mL/min [43.2-179.1], tumor size median 50.09 mm [8.14-305.24], soluble PD-L1 median 138.34 pg/mL [10.42-1478.7], anti-drug antibodies 0 for all subjects, height median 167.45 cm [143.5-198]. The underlying PK parameters were estimated by Baverel 2018 in a pooled oncology population (see the reference field)."
  )

  ini({
    # ---- Structural PK ---------------------------------------------------
    # All structural values are hardcoded literals in the ESM 3 $PK block of
    # a $SIMULATION ONLYSIM problem (the sole $THETA is "1 FIX", a dummy that
    # multiplies Q). Nothing here was estimated by de Vries 2025, so every
    # parameter is encoded as fixed(). CL, Q and Vmax are per DAY, per the
    # "; -- CL, Q and VMAX ARE PER DAY" comment in the control stream.
    lcl <- fixed(log(0.232)); label("Linear clearance at reference covariates (L/day)")                                # ESM 3 $PK: CL = 0.232 * (...)
    lvc <- fixed(log(3.51)); label("Central volume of distribution at reference covariates (L)")                       # ESM 3 $PK: V1 = 3.51 * (...)
    lvp <- fixed(log(3.45)); label("Peripheral volume of distribution at reference covariates (L)")                    # ESM 3 $PK: V2 = 3.45 * (...)
    lq <- fixed(log(0.476)); label("Intercompartmental clearance (L/day)")                                             # ESM 3 $PK: Q = 0.476 * THETA(1), THETA(1) = 1 FIX
    lvmax <- fixed(log(0.824)); label("Michaelis-Menten maximum elimination rate at reference sPD-L1 (mg/day)")        # ESM 3 $PK: VM = 0.824 * (...)
    lkm <- fixed(log(0.344)); label("Michaelis-Menten constant Km (mg/L)")                                             # ESM 3 $PK: KM = 0.344
    ld1 <- fixed(log(1 / 24)); label("Duration of the intravenous infusion D1 (day)")                                  # ESM 3 $PK: D1 = (1 / 24) ; 1 H INFUSION

    # ---- Covariate effects on linear clearance ---------------------------
    # The centred-linear terms are printed in ESM 3 as (1 - 0.035 * (ALB - 38))
    # and (1 + 0.00149 * (CRCL - 85.65)); both are encoded below in the
    # (1 + coefficient * (covariate - centre)) form, so the albumin slope
    # carries the minus sign.
    e_alb_cl <- fixed(-0.035); label("Slope of serum albumin on CL (per g/L above 38 g/L)")                            # ESM 3 $PK: (1 - 0.035 * (ALB - 38))
    e_crcl_cl <- fixed(0.00149); label("Slope of creatinine clearance on CL (per mL/min above 85.65 mL/min)")          # ESM 3 $PK: (1 + 0.00149 * (CRCL - 85.65))
    e_ecog_cl <- fixed(0.937); label("Multiplicative effect of ECOG performance status >= 1 on CL")                    # ESM 3 $PK: (0.937 ** FLAECOG)
    e_sexf_cl <- fixed(0.857); label("Multiplicative effect of female sex on CL")                                      # ESM 3 $PK: (0.857 ** FLASEX)
    e_tumsz_cl <- fixed(0.00178); label("Slope of baseline tumor size on CL (per mm above 74.8 mm)")                   # ESM 3 $PK: (1 + 0.00178 * (TUMORSIZE - 74.8))
    e_wt_cl <- fixed(0.389); label("Allometric exponent of body weight on CL (unitless)")                              # ESM 3 $PK: ((WT / 69.8) ** 0.389)

    # ---- Covariate effects on Vmax and the volumes -----------------------
    e_spdl1_vmax <- fixed(0.00336); label("Slope of soluble PD-L1 on Vmax (per pg/mL above 124.8 pg/mL)")              # ESM 3 $PK: (1 + 0.00336 * (SPDL1 - 124.8))
    e_sexf_vc <- fixed(0.835); label("Multiplicative effect of female sex on Vc")                                      # ESM 3 $PK: V1 = 3.51 * (0.835 ** FLASEX) * (...)
    e_wt_vc <- fixed(0.406); label("Allometric exponent of body weight on Vc (unitless)")                              # ESM 3 $PK: ((WT / 69.8) ** 0.406)
    e_sexf_vp <- fixed(0.795); label("Multiplicative effect of female sex on Vp")                                      # ESM 3 $PK: V2 = 3.45 * (0.795 ** FLASEX) * (...)

    # ---- Between-subject variability -------------------------------------
    # ESM 3 "$OMEGA BLOCK(3) CORRELATION" with a FIX flag on the block:
    #   diagonal   variances 0.0729 (CL), 0.0437 (V1), 0.113 (V2)
    #   off-diagonal CORRELATIONS 0.279 (CL-V1), 0 (CL-V2), 0.560 (V1-V2)
    # Converted to covariances below as corr * sd_i * sd_j:
    #   cov(CL, V1) = 0.279 * sqrt(0.0729) * sqrt(0.0437) = 0.01574739
    #   cov(CL, V2) = 0
    #   cov(V1, V2) = 0.560 * sqrt(0.0437) * sqrt(0.113)  = 0.03935210
    # Approximate log-scale CVs: 27.5 percent (CL), 21.1 percent (Vc),
    # 34.6 percent (Vp).
    etalcl + etalvc + etalvp ~ fixed(c(
      0.0729,
      0.01574739, 0.0437,
      0, 0.03935210, 0.113
    ))                                                                                                                 # ESM 3 $OMEGA BLOCK(3) CORRELATION, entries 5-7

    # ---- Residual unexplained variability --------------------------------
    # ESM 3 $ERROR: Y = F + (F * ERR(1)) + ERR(2), i.e. combined proportional
    # plus additive on the linear concentration scale. $SIGMA holds VARIANCES,
    # so the standard deviations below are their square roots:
    #   propSd = sqrt(0.0454) = 0.2130728 (fraction)
    #   addSd  = sqrt(0.09)   = 0.3 mg/L
    propSd <- fixed(0.2130728); label("Proportional residual error (fraction)")                                        # ESM 3 $SIGMA 0.0454 ; PROP ERR
    addSd <- fixed(0.3); label("Additive residual error (mg/L)")                                                       # ESM 3 $SIGMA 0.09 ; ADD ERR
  })

  model({
    # ---- Individual parameters -------------------------------------------
    # Covariate forms transcribed verbatim from the ESM 3 $PK block. The
    # centring constants 38 g/L, 85.65 mL/min, 74.8 mm, 124.8 pg/mL and the
    # allometric reference 69.8 kg are the Baverel 2018 reference-patient
    # values printed inside those expressions.
    cl <- exp(lcl + etalcl) *
      (1 + e_alb_cl * (ALB - 38)) *
      (1 + e_crcl_cl * (CRCL - 85.65)) *
      e_ecog_cl^ECOG_GE1 *
      e_sexf_cl^SEXF *
      (1 + e_tumsz_cl * (TUMSZ - 74.8)) *
      (WT / 69.8)^e_wt_cl
    vmax <- exp(lvmax) * (1 + e_spdl1_vmax * (SPDL1 - 124.8))
    vc <- exp(lvc + etalvc) * e_sexf_vc^SEXF * (WT / 69.8)^e_wt_vc
    vp <- exp(lvp + etalvp) * e_sexf_vp^SEXF
    q <- exp(lq)
    km <- exp(lkm)
    d1 <- exp(ld1)

    # ---- Micro-constants (ESM 3 $PK: K10 = CL/V1, K12 = Q/V1, K21 = Q/V2) --
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- Plasma concentration (ESM 3 $PK sets S1 = V1, so F = A(1)/V1) -----
    Cc <- central / vc

    # ---- ODE system -------------------------------------------------------
    # ESM 3 $DES, verbatim:
    #   DADT(1) = -K10*A(1) - K12*A(1) - ((VM*C1)/(KM + C1)) + K21*A(2)
    #   DADT(2) = -K21*A(2) + K12*A(1)
    # Linear (FcRn-mediated) elimination and a saturable Michaelis-Menten arm
    # act in parallel from the central compartment. The ESM 3 $MODEL block
    # also declares AUCSS / AUC1 / AUCTOTAL accumulator compartments whose
    # DADT are C1 gated on fixed calendar windows (T >= 252 d, T < 28 d);
    # those are reporting bookkeeping for the paper's own regimens, carry no
    # PK information and would hard-code that paper's timeline into a reusable
    # library model, so they are computed post hoc in the vignette instead.
    d/dt(central) <- -kel * central - k12 * central -
      (vmax * Cc) / (km + Cc) + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # ---- Zero-order intravenous infusion ----------------------------------
    # Dose rows must carry rate = -2 for rxode2 to use the modelled duration.
    dur(central) <- d1

    Cc ~ add(addSd) + prop(propSd)
  })
}
