Lacy_2018_cabozantinib <- function() {
  description <- "Two-compartment population PK model for oral cabozantinib (tyrosine kinase inhibitor) in healthy volunteers and patients with renal cell carcinoma, castration-resistant prostate cancer, medullary thyroid carcinoma, glioblastoma multiforme, or other advanced malignancies (Lacy 2018, n=1534 across 9 clinical studies). Absorption is described by parallel dual processes: a fraction F1 enters depot1 via first-order absorption rate Ka with absorption lag time ALAG1, and the remaining (1-F1) enters depot2 via zero-order infusion over duration D2. Capsule (vs tablet, reference) formulation reduces both Ka and overall bioavailability; Ka also scales with dose via a power function (DOSE/60 mg)^0.677. Two-compartment disposition (central + peripheral1) with first-order elimination from central. Covariates on CL/F and Vc/F are baseline age (power on median 64 y), body weight (power on median 81 kg), female sex, race (Black/Asian/Other vs White reference), and tumor type (RCC/CRPC/MTC/GB/Other vs HV reference); MTC cancer type drives an approximately 93% higher CL/F."
  reference <- paste(
    "Lacy S, Yang B, Nielsen J, Miles D, Nguyen L, Hutmacher M.",
    "A population pharmacokinetic model of cabozantinib in healthy volunteers",
    "and patients with various cancer types.",
    "Cancer Chemother Pharmacol. 2018;81(6):1071-1082.",
    "doi:10.1007/s00280-018-3581-0.",
    sep = " "
  )
  vignette <- "Lacy_2018_cabozantinib"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    AGE = list(
      description        = "Subject age at baseline",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline. Power effect on CL/F (exponent -0.162) and Vc/F (exponent -0.012). Reference 64 years is the overall-cohort median across the nine pooled studies (Lacy 2018 Table 2). Paper Methods: 'The approximate median value was used for xREF.'",
      source_name        = "AGE"
    ),
    WT = list(
      description        = "Body weight at baseline",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline. Power effect on CL/F (exponent -0.028, near-null) and Vc/F (exponent +1.019, near-linear scaling). Reference 81 kg is the overall-cohort median across the nine pooled studies (Lacy 2018 Table 2). Paper Methods: 'The approximate median value was used for xREF.'",
      source_name        = "WT"
    ),
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male, the typical-value reference)",
      notes              = "Time-fixed. Multiplicative fractional effect on CL/F (-0.230, i.e. 23% lower in females) and Vc/F (+0.11). Lacy 2018 Table 3.",
      source_name        = "SEXF"
    ),
    RACE_BLACK = list(
      description        = "Black / African American race indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Black; White is the typical-value reference category when paired with RACE_ASIAN = 0 and RACE_OTHER = 0)",
      notes              = "Time-fixed. Multiplicative fractional effect on CL/F (+0.301, i.e. 30% higher in Black subjects) and Vc/F (-0.022). Lacy 2018 Table 3. Reference = White (1283 / 1534 = 83.6% of cohort).",
      source_name        = "RACE"
    ),
    RACE_ASIAN = list(
      description        = "Asian race indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Asian; White is the typical-value reference when paired with RACE_BLACK = 0 and RACE_OTHER = 0)",
      notes              = "Time-fixed. Multiplicative fractional effect on CL/F (-0.078) and Vc/F (+0.05). Lacy 2018 Table 3.",
      source_name        = "RACE"
    ),
    RACE_OTHER = list(
      description        = "Race category 'Other' indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Other; White is the typical-value reference when paired with RACE_BLACK = 0 and RACE_ASIAN = 0)",
      notes              = "Time-fixed. Multiplicative fractional effect on CL/F (-0.007, near-null) and Vc/F (-0.059). Lacy 2018 Table 3. Race 'Other' includes mixed and unspecified categories pooled across the nine studies (39 / 1534 = 2.5% of cohort).",
      source_name        = "RACE"
    ),
    TUMTP_RCC = list(
      description        = "Renal cell carcinoma indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-RCC; healthy volunteer is the typical-value reference when paired with all other TUMTP_* indicators = 0)",
      notes              = "Time-fixed. Multiplicative fractional effect on CL/F (-0.129, 12.9% lower than HV) and Vc/F (-0.63, 63% lower than HV). Lacy 2018 Table 3. Reference = healthy volunteer (140 / 1534 = 9.1% of cohort).",
      source_name        = "POP"
    ),
    TUMTP_HRPC = list(
      description        = "Castration-resistant prostate cancer indicator (paper writes CRPC; canonical column TUMTP_HRPC covers both HRPC and CRPC wordings)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-CRPC; healthy volunteer is the typical-value reference)",
      notes              = "Time-fixed. Multiplicative fractional effect on CL/F (-0.009, near-null) and Vc/F (-0.241). Lacy 2018 Table 3. Reference = healthy volunteer. CRPC patients accounted for 823 / 1534 = 53.7% of the pooled cohort across Studies 203, 306, and 307.",
      source_name        = "POP"
    ),
    TUMTP_MTC = list(
      description        = "Medullary thyroid carcinoma indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-MTC; healthy volunteer is the typical-value reference)",
      notes              = "Time-fixed. Multiplicative fractional effect on CL/F (+0.928, 92.8% higher than HV; load-bearing finding of Lacy 2018) and Vc/F (-0.07). Reference = healthy volunteer. MTC accounted for 210 / 1534 = 13.7% of cohort, all from Study 301 dosed at 140 mg/day capsule. The MTC effect on CL/F is reported as time-dependent (negligible on day 1, fully expressed by day 29).",
      source_name        = "POP"
    ),
    TUMTP_GLIO = list(
      description        = "Glioblastoma multiforme indicator (canonical TUMTP_GLIO covers glioma of any grade including GB)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-glioma; healthy volunteer is the typical-value reference)",
      notes              = "Time-fixed. Multiplicative fractional effect on CL/F (+0.216) and Vc/F (-0.569). Lacy 2018 Table 3. Reference = healthy volunteer. GB patients accounted for 39 / 1534 = 2.5% of the cohort (Study 201).",
      source_name        = "POP"
    ),
    TUMTP_OTH = list(
      description        = "Heterogeneous 'other malignancy' pool indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-other; healthy volunteer is the typical-value reference)",
      notes              = "Time-fixed. Multiplicative fractional effect on CL/F (+0.178) and Vc/F (-0.186). Lacy 2018 Table 3. Reference = healthy volunteer. 'Other' = mixed advanced malignancies enrolled in the FIH Study 001 (n = 40 / 1534 = 2.6% of cohort); per-subject tumor composition is not enumerated by Lacy 2018.",
      source_name        = "POP"
    ),
    FORM_CAPSULE = list(
      description        = "Capsule formulation indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (tablet; the typical-value reference)",
      notes              = "Per-dose-occasion indicator (Study 010 is a tablet-vs-capsule crossover). Multiplicative fractional effect on Ka (-0.579, i.e. 57.9% lower absorption rate for capsule) and on overall bioavailability F (-0.144, i.e. 14.4% lower exposure for capsule). Reference = tablet (Cabometyx 60 mg approved for RCC / CRPC). Comparator capsule = Cometriq 140 mg approved for MTC. Aligns with the cross-study capsule-tablet bioequivalence finding (Study XL184-010).",
      source_name        = "FORM"
    ),
    DOSE = list(
      description        = "Administered cabozantinib dose level",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Per-dose-occasion. Power covariate on first-order absorption rate constant Ka: Ka(DOSE) = Ka_ref * (DOSE/DOSE_REF)^0.677. Reference DOSE_REF = 60 mg (the standard tablet daily dose for non-MTC indications and the reference dose used in the companion Lacy 2018 exposure-response paper). The DOSE_REF value is NOT stated in the paper; the 60 mg choice was operator-approved per the sidecar resolution Q2-A and is documented in vignette Errata.",
      source_name        = "DOSE"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 1534L,
    n_observations   = 8072L,
    n_studies        = 9L,
    age_range        = "18-87 years",
    age_median       = "64 years (pooled across the nine studies)",
    weight_range     = "30.4-190.7 kg",
    weight_median    = "approximately 81 kg (pooled across the nine studies)",
    sex_female_pct   = 14.4,
    race_ethnicity   = "White 83.6%, Black 2.7%, Asian 3.0%, Other 2.5%, missing/not reported 8.1% (Lacy 2018 Table 2)",
    disease_state    = "Pooled population: healthy volunteers (9.1%); castration-resistant prostate cancer (53.7%); renal cell carcinoma (18.4%); medullary thyroid cancer (13.7%); glioblastoma multiforme (2.5%); other advanced malignancies (2.6%, Study 001 mixed cohort).",
    dose_range       = "20-200 mg/day oral cabozantinib (capsule Cometriq 140 mg/day for MTC; tablet Cabometyx 60 mg/day for RCC / CRPC; phase I HV studies 20, 40, 60, 140 mg single doses).",
    regions          = "Multinational (Lacy 2018 Table 1 enumerates the nine studies; specific regional breakdown not reported)",
    formulations     = "Capsule (Cometriq, 42.2% of cohort) and tablet (Cabometyx, 62.5%); a subset of subjects in the bioequivalence Study 010 received both formulations in a crossover design.",
    studies          = c("XL184-001 (FIH, mixed malignancies, 140 or 200 mg, n=40)",
                         "XL184-010 (HV BE crossover capsule vs tablet, 140 mg, n=77)",
                         "XL184-020 (HV tablet 20 / 40 / 60 mg, n=63)",
                         "XL184-201 (GB, 140 mg QD, n=39)",
                         "XL184-203 (CRPC RDT / NRE, 40 or 100 mg QD, n=284)",
                         "XL184-301 (MTC, 140 mg QD, n=210)",
                         "XL184-306 (CRPC, 60 mg QD, n=41)",
                         "XL184-307 (CRPC, 60 mg QD, n=498)",
                         "XL184-308 (RCC, 60 mg QD, n=282)"),
    notes            = "Baseline demographics from Lacy 2018 Table 2. BLQ rate < 1% (excluded from analysis per Methods). 79% of the 210 MTC patients in Study 301 dose-reduced from the starting 140 mg. The MTC subjects had ~93% higher CL/F than HV at steady state but not at day 1, suggesting time-varying clearance in this subgroup; the FM model captures the steady-state difference."
  )

  ini({
    # ---- Structural population parameters (Lacy 2018 Table 3, "Full model (FM)" column) ----
    # Reference covariate set: tablet formulation, 60 mg dose, healthy volunteer
    # (HV reference), male, White, age 64 y, weight 81 kg (overall-cohort
    # medians from Table 2). Apparent (oral) parameters CL/F, Vc/F, Q/F, Vp/F
    # are reported in L/h and L; Ka in 1/h; ALAG1 and D2 in h.
    lka     <- log(0.979)  ; label("First-order absorption rate constant at the 60 mg tablet reference (Ka, 1/h)")               # Lacy 2018 Table 3 FM Ka = 0.979 (90% CI 0.679, 1.411)
    ld2     <- log(2.4)    ; label("Zero-order absorption duration for the parallel depot2 process (D2, h)")                     # Lacy 2018 Table 3 FM D2 = 2.4 (90% CI 2.01, 2.866)
    lcl     <- log(2.478)  ; label("Apparent oral clearance at the reference covariate set (CL/F, L/h)")                         # Lacy 2018 Table 3 FM CL/F = 2.478 (90% CI 2.257, 2.721)
    lvc     <- log(187.0)  ; label("Apparent central volume of distribution at the reference covariate set (Vc/F, L)")           # Lacy 2018 Table 3 FM Vc/F = 187.0 (90% CI 156.3, 223.9)
    lq      <- log(31.213) ; label("Apparent inter-compartmental clearance (Q/F, L/h)")                                          # Lacy 2018 Table 3 FM Q/F = 31.213 (90% CI 28.732, 33.92)
    lvp     <- log(195.1)  ; label("Apparent peripheral volume of distribution (Vp/F, L)")                                       # Lacy 2018 Table 3 FM Vp/F = 195.1 (90% CI 183.3, 207.9)
    lalag1  <- log(0.784)  ; label("Absorption lag time for the parallel depot1 first-order process (ALAG1, h)")                 # Lacy 2018 Table 3 FM ALAG1 = 0.784 (90% CI 0.757, 0.812)

    # F1 = fraction of the dose absorbed via the first-order depot1; the
    # remaining (1 - F1) is absorbed via the zero-order depot2 process. Paper
    # footnote a: "Anti-logit transformation was used to obtain F1", so the
    # estimated THETA lives on the logit scale and the transformed value 0.854
    # is the back-transformed proportion. logit(0.854) = log(0.854 / 0.146)
    # = 1.7665.
    logitf1 <- 1.7665      ; label("Logit of the fraction absorbed via the first-order depot1 process (F1 = expit(logitf1))")    # Lacy 2018 Table 3 FM F1 (transformed) = 0.854 (90% CI 0.819, 0.884); logit(0.854) = 1.7665

    # Bioavailability anchor. The reference (tablet) overall F is fixed at 1
    # because Lacy 2018 reports the capsule-vs-tablet F covariate only
    # (no separate F estimate at the tablet reference).
    lfdepot <- fixed(log(1.0)) ; label("Overall oral bioavailability at the tablet reference (F, fixed at 1)")                   # Lacy 2018 Methods; tablet is the reference formulation

    # ---- Dose-dependent Ka power exponent ----
    # Paper Results: "The first-order absorption process including ... a dose-
    # dependent effect on the absorption rate constant (Ka) was described using
    # a power model." Ka(DOSE) = Ka_ref * (DOSE / DOSE_REF)^exp_dose_ka.
    # DOSE_REF is NOT stated in the paper; the 60 mg reference was chosen per
    # sidecar Q2-A and documented in vignette Errata.
    e_dose_ka <- 0.677     ; label("Dose power exponent on Ka (unitless)")                                                       # Lacy 2018 Table 3 FM dose-dependent Ka = 0.677 (90% CI 0.268, 1.085)

    # ---- Categorical covariate effects (multiplicative fractional change) ----
    # Paper footnote b: "For categorical covariates ... transformed estimates
    # correspond to fractional change from the reference level." Encoded as
    # TV_with_cov = TV_ref * (1 + e_<cov>_<param> * IND), giving the published
    # fractional change at IND = 1.

    # Formulation effects (capsule vs tablet)
    e_form_capsule_ka <- -0.579 ; label("Capsule (vs tablet) fractional change on Ka (unitless)")                                # Lacy 2018 Table 3 FM Capsule on Ka = -0.579 (90% CI -0.783, -0.183)
    e_form_capsule_f  <- -0.144 ; label("Capsule (vs tablet) fractional change on overall oral bioavailability F (unitless)")    # Lacy 2018 Table 3 FM Capsule on overall relative oral availability = -0.144 (90% CI -0.162, -0.126)

    # Sex effects (female vs male reference)
    e_sexf_cl <- -0.230 ; label("Female (vs male) fractional change on CL/F (unitless)")                                         # Lacy 2018 Table 3 FM Female on CL/F = -0.230 (90% CI -0.286, -0.17)
    e_sexf_vc <-  0.11  ; label("Female (vs male) fractional change on Vc/F (unitless)")                                         # Lacy 2018 Table 3 FM Female on Vc/F = 0.11 (90% CI -0.033, 0.276)

    # Race effects (Black / Asian / Other vs White reference)
    e_race_black_cl <-  0.301 ; label("Race Black (vs White) fractional change on CL/F (unitless)")                              # Lacy 2018 Table 3 FM Race (Black) on CL/F = 0.301
    e_race_black_vc <- -0.022 ; label("Race Black (vs White) fractional change on Vc/F (unitless)")                              # Lacy 2018 Table 3 FM Race (Black) on Vc/F = -0.022
    e_race_asian_cl <- -0.078 ; label("Race Asian (vs White) fractional change on CL/F (unitless)")                              # Lacy 2018 Table 3 FM Race (Asian) on CL/F = -0.078
    e_race_asian_vc <-  0.05  ; label("Race Asian (vs White) fractional change on Vc/F (unitless)")                              # Lacy 2018 Table 3 FM Race (Asian) on Vc/F = 0.05
    e_race_other_cl <- -0.007 ; label("Race Other (vs White) fractional change on CL/F (unitless)")                              # Lacy 2018 Table 3 FM Race (Other) on CL/F = -0.007
    e_race_other_vc <- -0.059 ; label("Race Other (vs White) fractional change on Vc/F (unitless)")                              # Lacy 2018 Table 3 FM Race (Other) on Vc/F = -0.059

    # Cancer-type effects (RCC / CRPC / MTC / GB / Other vs HV reference)
    e_tumtp_rcc_cl  <- -0.129 ; label("RCC (vs HV) fractional change on CL/F (unitless)")                                        # Lacy 2018 Table 3 FM RCC on CL/F = -0.129
    e_tumtp_rcc_vc  <- -0.63  ; label("RCC (vs HV) fractional change on Vc/F (unitless)")                                        # Lacy 2018 Table 3 FM RCC on Vc/F = -0.63
    e_tumtp_hrpc_cl <- -0.009 ; label("CRPC (vs HV) fractional change on CL/F (unitless)")                                       # Lacy 2018 Table 3 FM CRPC on CL/F = -0.009
    e_tumtp_hrpc_vc <- -0.241 ; label("CRPC (vs HV) fractional change on Vc/F (unitless)")                                       # Lacy 2018 Table 3 FM CRPC on Vc/F = -0.241
    e_tumtp_mtc_cl  <-  0.928 ; label("MTC (vs HV) fractional change on CL/F (unitless)")                                        # Lacy 2018 Table 3 FM MTC on CL/F = 0.928
    e_tumtp_mtc_vc  <- -0.07  ; label("MTC (vs HV) fractional change on Vc/F (unitless)")                                        # Lacy 2018 Table 3 FM MTC on Vc/F = -0.07
    e_tumtp_glio_cl <-  0.216 ; label("GB (vs HV) fractional change on CL/F (unitless)")                                         # Lacy 2018 Table 3 FM GB on CL/F = 0.216
    e_tumtp_glio_vc <- -0.569 ; label("GB (vs HV) fractional change on Vc/F (unitless)")                                         # Lacy 2018 Table 3 FM GB on Vc/F = -0.569
    e_tumtp_oth_cl  <-  0.178 ; label("Other malignancies (vs HV) fractional change on CL/F (unitless)")                         # Lacy 2018 Table 3 FM Other malignancies on CL/F = 0.178
    e_tumtp_oth_vc  <- -0.186 ; label("Other malignancies (vs HV) fractional change on Vc/F (unitless)")                         # Lacy 2018 Table 3 FM Other malignancies on Vc/F = -0.186

    # ---- Continuous covariate effects (power exponents, untransformed) ----
    e_age_cl <- -0.162 ; label("Age power exponent on CL/F (AGE / 64 y; unitless)")                                              # Lacy 2018 Table 3 FM Age on CL/F = -0.162
    e_age_vc <- -0.012 ; label("Age power exponent on Vc/F (AGE / 64 y; unitless)")                                              # Lacy 2018 Table 3 FM Age on Vc/F = -0.012
    e_wt_cl  <- -0.028 ; label("Body-weight power exponent on CL/F (WT / 81 kg; unitless)")                                      # Lacy 2018 Table 3 FM Weight on CL/F = -0.028
    e_wt_vc  <-  1.019 ; label("Body-weight power exponent on Vc/F (WT / 81 kg; unitless)")                                      # Lacy 2018 Table 3 FM Weight on Vc/F = 1.019

    # ---- Inter-individual variability ----
    # Lacy 2018 Table 3 footnote d reports untransformed variances:
    #   sigma^2 = 0.118; omega^2_Ka = 2.063; omega^2_CL/F = 0.202;
    #   omega^2_Vc/F = 0.233; omega^2_F1 = 0.466 (on the logit scale);
    #   omega^2_CL/F:Vc/F = 2.475.
    # The published off-diagonal value 2.475 violates Cauchy-Schwarz given the
    # marginal variances (implied correlation 11.4); per sidecar Q1-A the
    # off-diagonal is dropped and Omega is encoded as diagonal. See vignette
    # Errata.
    etalka      ~ 2.063   # Lacy 2018 Table 3 footnote d omega^2_Ka = 2.063
    etalcl      ~ 0.202   # Lacy 2018 Table 3 footnote d omega^2_CL/F = 0.202; off-diagonal with etalvc dropped (sidecar Q1-A; see vignette Errata)
    etalvc      ~ 0.233   # Lacy 2018 Table 3 footnote d omega^2_Vc/F = 0.233; off-diagonal with etalcl dropped (sidecar Q1-A; see vignette Errata)
    etalogitf1  ~ 0.466   # Lacy 2018 Table 3 footnote d omega^2_F1 = 0.466 (on the logit scale, matching the F1 estimation scale per footnote a)

    # ---- Residual error ----
    # Paper Methods: "Residual variability ... was initially modeled using the
    # log-transformed additive error model". An additive residual on
    # log-transformed concentration in NONMEM (LTBS) corresponds to a
    # proportional residual in linear-space nlmixr2. sqrt(0.118) = 0.3435.
    propSd <- 0.3435 ; label("Proportional residual error (fraction; sqrt(sigma^2 = 0.118) from the LTBS additive error)")       # Lacy 2018 Table 3 footnote d sigma^2 = 0.118
  })

  model({
    # ---- Reference covariate values (Lacy 2018 Methods: 'The approximate
    # median value was used for xREF.'  Overall-cohort medians from Table 2.) ----
    ref_age  <- 64    # years
    ref_wt   <- 81    # kg
    ref_dose <- 60    # mg (sidecar Q2-A; see vignette Errata)

    # ---- Categorical covariate multipliers (multiplicative fractional change) ----
    # TV_with_cov = TV_ref * (1 + e_<cov>_<param> * IND); reference category
    # gives multiplier 1.
    cl_sex   <- 1 + e_sexf_cl * SEXF
    vc_sex   <- 1 + e_sexf_vc * SEXF

    cl_race  <- 1 +
                e_race_black_cl * RACE_BLACK +
                e_race_asian_cl * RACE_ASIAN +
                e_race_other_cl * RACE_OTHER
    vc_race  <- 1 +
                e_race_black_vc * RACE_BLACK +
                e_race_asian_vc * RACE_ASIAN +
                e_race_other_vc * RACE_OTHER

    cl_tumtp <- 1 +
                e_tumtp_rcc_cl  * TUMTP_RCC  +
                e_tumtp_hrpc_cl * TUMTP_HRPC +
                e_tumtp_mtc_cl  * TUMTP_MTC  +
                e_tumtp_glio_cl * TUMTP_GLIO +
                e_tumtp_oth_cl  * TUMTP_OTH
    vc_tumtp <- 1 +
                e_tumtp_rcc_vc  * TUMTP_RCC  +
                e_tumtp_hrpc_vc * TUMTP_HRPC +
                e_tumtp_mtc_vc  * TUMTP_MTC  +
                e_tumtp_glio_vc * TUMTP_GLIO +
                e_tumtp_oth_vc  * TUMTP_OTH

    ka_form  <- 1 + e_form_capsule_ka * FORM_CAPSULE
    f_form   <- 1 + e_form_capsule_f  * FORM_CAPSULE

    # ---- Individual PK parameters ----
    ka <- exp(lka + etalka) *
          (DOSE / ref_dose)^e_dose_ka *
          ka_form
    cl <- exp(lcl + etalcl) *
          (AGE / ref_age)^e_age_cl *
          (WT  / ref_wt )^e_wt_cl  *
          cl_sex * cl_race * cl_tumtp
    vc <- exp(lvc + etalvc) *
          (AGE / ref_age)^e_age_vc *
          (WT  / ref_wt )^e_wt_vc  *
          vc_sex * vc_race * vc_tumtp
    q  <- exp(lq)
    vp <- exp(lvp)

    d2     <- exp(ld2)
    alag1  <- exp(lalag1)
    f1     <- expit(logitf1 + etalogitf1)
    fdepot <- exp(lfdepot) * f_form

    # ---- Micro-constants for the explicit ODE system ----
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ---- ODE system: parallel dual absorption + 2-compartment disposition ----
    # depot1 receives a fraction F1 of the dose with first-order absorption
    # rate Ka and lag ALAG1; the remaining (1 - F1) of the dose is
    # zero-order-infused directly into central over duration D2 (no holding
    # depot). Bioavailability F (capsule vs tablet) multiplies both fractions.
    # Dosing events must target BOTH compartments (cmt = depot1, cmt = central)
    # with the same total amount; the f() multipliers split the dose into the
    # F1 and (1 - F1) fractions and the dur(central) on the central-targeted
    # row produces the zero-order infusion. See the vignette make_cohort for
    # the event-table construction pattern.
    d/dt(depot1)      <- -ka * depot1
    d/dt(central)     <-  ka * depot1 - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    f(depot1)    <- fdepot * f1
    f(central)   <- fdepot * (1 - f1)
    alag(depot1) <- alag1
    dur(central) <- d2

    # ---- Plasma concentration ----
    # central is in mg, vc in L -> mg/L; multiply by 1000 to report in ng/mL
    # (matches the Lacy 2018 bioanalytical LLOQ unit, 0.5 ng/mL).
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
