Lu_2019_polatuzumab <- function() {
  description <- "Integrated two-analyte population PK model of polatuzumab vedotin (anti-CD79b vc-MMAE antibody-drug conjugate) in adults with non-Hodgkin lymphoma (Lu 2019). The antibody-conjugated MMAE (acMMAE) is described by a two-compartment model with three parallel elimination pathways from the central compartment: a slowly-time-decaying nonspecific linear clearance (CL_NS, sigmoidal Hill decline with cycle), a rapidly-decaying linear clearance (CL_t, mono-exponential decline), and a saturable Michaelis-Menten clearance (CL_MM). All three acMMAE pathways feed unconjugated MMAE formation in the central MMAE compartment with relative conversion fractions FRAC_NS, FRAC_NS x FRAC_CLT, and FRAC_NS x FRAC_MM, modulated by a time-dependent multiplier (1 + FRAC_T x exp(-alpha x t)) on FRAC_NS that captures the cycle-over-cycle decline in MMAE formation. Unconjugated MMAE is described by an apparent two-compartment model with parallel linear (CL_MMAE) and Michaelis-Menten (Vmax_MMAE / KSS) elimination from its central compartment. Modeled in MMAE-equivalent micrograms (pola dose in ug/kg x weight in kg x 3.65 x 718 / 145001 -> MMAE-equivalent ug administered to the acMMAE central compartment), with concentrations in ng/mL = ug/L."
  reference <- "Lu D, Lu T, Gibiansky L, Li X, Li C, Agarwal P, Shemesh CS, Shi R, Dere RC, Hirata J, Miles D, Chanu P, Girish S, Jin JY. Integrated Two-Analyte Population Pharmacokinetic Model of Polatuzumab Vedotin in Patients With Non-Hodgkin Lymphoma. CPT Pharmacometrics Syst Pharmacol. 2020;9(1):48-59. doi:10.1002/psp4.12482. PMID 31749251."
  vignette <- "Lu_2019_polatuzumab"
  units <- list(
    time          = "hour",
    dosing        = "ug",
    concentration = "ng/mL"
  )

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline value. Power effects on acMMAE CL_SS (exponent 0.73), on Vc/Vp/Q (shared exponent 0.50), and on FRAC_NS (exponent -0.467); reference 75 kg per Lu 2019 NONMEM control stream (BWT/75 normalization).",
      source_name        = "BWT"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Time-fixed. Lu 2019 NONMEM defines SEX = SEXN - 1 (SEXN = 1 female, 2 male) so the source 'SEX' indicator equals 1 for males. The canonical SEXF reverses the value coding (SEXF = 1 - source_SEX); effect-coefficient ratios therefore invert (e.g., paper reports V1_male/V1_female = 1.20, stored as e_sexf_vc = 1/1.20 = 0.8333 applied as e_sexf_vc^SEXF). See Assumptions section of the validation vignette for the full sign-and-reference-category derivation.",
      source_name        = "SEXN"
    ),
    LINE_1L = list(
      description        = "First-line-therapy indicator: 1 = treatment-naive (previously untreated), 0 = relapsed/refractory (>= second line)",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Time-fixed. Lu 2019 NONMEM defines NAIVE = 1 if RRFN == 0 (treatment-naive); the canonical LINE_1L matches this coding directly (no value flip). Treatment-naive status enters as multiplicative effects on Vc, kdes (the rate constant of CL_TIME decay), CL_TIME (initial linear time-decaying clearance), and FRAC_NS.",
      source_name        = "RRFN"
    ),
    RACE_ASIAN = list(
      description        = "Asian race indicator, 1 = Asian, 0 = non-Asian",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Time-fixed. Lu 2019 NONMEM defines ASIAN = 1 if RACEN == 1; the canonical RACE_ASIAN matches this coding (no value flip). Multiplicative effect on Vc only (e_asian_vc = 0.929, applied as e_asian_vc^RACE_ASIAN).",
      source_name        = "RACEN"
    ),
    ALB = list(
      description        = "Baseline serum albumin concentration",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline value (SI units). Power effects on CL_SS (exponent -0.247) and FRAC_NS (exponent -0.613); reference 35 g/L per Lu 2019 NONMEM control stream (BALBUM/35 normalization).",
      source_name        = "BALBUM"
    ),
    BLBCELL = list(
      description        = "Baseline CD19+ B cell count",
      units              = "cells/uL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline value. Lu 2019 reports B-cell count in 10^6 cells/L which equals cells/uL by unit conversion (10^6 / L = 1 / uL). Two distinct effects: (1) on CL_SS, applied as max(1, BLBCELL)^0.0212 (i.e., a power effect with the input floored at 1 cell/uL); (2) on CL_TIME, applied as max(1, BLBCELL/121)^0.578 (a power effect with the threshold of 121 cells/uL below which the multiplier is 1).",
      source_name        = "BBCC"
    ),
    TUMSZ = list(
      description        = "Baseline tumor sum of the products of perpendicular dimensions (SPD) for the non-Hodgkin lymphoma cohort",
      units              = "mm^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline value. The canonical TUMSZ entry pools SPD (NHL/cHL, mm^2) and RECIST sum-of-diameters (solid tumors, mm) onto one column; for Lu 2019 the unit is mm^2 (SPD) per the source paper's NHL methodology. Two distinct effects: (1) on CL_SS, applied as 1 + 0.0521 x (TUMSZ/5000 - 1) (linear effect normalized to 5000 mm^2 reference, equivalent to NONMEM 1 + theta31*(BTMBD/5000-1)); (2) on CL_TIME, applied as TUMSZ / (1150 + TUMSZ) (Michaelis-Menten / saturable form with 50% effect at 1150 mm^2). Sensitivity-analysis range observed: 355 mm^2 (5th percentile) to 20000 mm^2 (95th percentile).",
      source_name        = "BTMBD"
    ),
    ECOG_GE1 = list(
      description        = "Eastern Cooperative Oncology Group (ECOG) performance status indicator, 1 = ECOG >= 1, 0 = ECOG = 0",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Time-fixed. Lu 2019 NONMEM defines ECOG0 = 1 if BECOG == 0; the canonical ECOG_GE1 reverses the value coding (ECOG_GE1 = 1 - source_ECOG0) and the effect-coefficient ratio therefore inverts (paper reports FRAC_NS_ECOG=0 / FRAC_NS_ECOG>=1 = 0.905, stored as e_ecog_ge1_frac_mmae = 1/0.905 = 1.1050 applied as e_ecog_ge1_frac_mmae^ECOG_GE1). Effect on FRAC_NS only.",
      source_name        = "BECOG"
    ),
    HEPIMP = list(
      description        = "Baseline hepatic-impairment indicator per the National Cancer Institute Organ Dysfunction Working Group (NCI ODWG) classification, 1 = mild or worse hepatic impairment (NCI ODWG group >= mild), 0 = normal hepatic function",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Time-fixed. Lu 2019 NONMEM defines HEPA = 1 if BHPTGRPN > 1.5 (i.e., NCI ODWG group >= 2 = mild or worse) AND BHPTGRPN != 9999 (the missing-value sentinel). Multiplicative effect on FRAC_NS only (e_hepimp_frac_mmae = 1.19 applied as e_hepimp_frac_mmae^HEPIMP). NCI ODWG group 1 = normal (reference); group 2 = mild; group 3 = moderate; group 4 = severe (Ramalingam et al., J Clin Oncol 2010;28:4507).",
      source_name        = "BHPTGRPN"
    ),
    COMBO_RG = list(
      description        = "Anti-CD20 combination-therapy indicator, 1 = polatuzumab vedotin co-administered with rituximab OR obinutuzumab, 0 = single-agent polatuzumab vedotin",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Time-fixed. Lu 2019 NONMEM defines RTX = 1 if COMBO == 1, GA101 = 1 if COMBO == 2, and applies effects as theta^(RTX+GA101); since RTX and GA101 are mutually exclusive, RTX+GA101 takes values 0 or 1 and the effect collapses to a single anti-CD20-combination indicator. Multiplicative effects on CL_SS (e_combo_rg_cl = 0.844, lower CL_SS on combo), kdes (e_combo_rg_kdes = 0.932), and FRAC_NS (e_combo_rg_frac_mmae = 0.709). The paper distinguishes rituximab and obinutuzumab combos in Table S1 but the final model fits a single combined effect (i.e., does not detect a meaningful difference between the two).",
      source_name        = "COMBO"
    )
  )

  population <- list(
    n_subjects     = 460L,
    n_studies      = 4L,
    n_observations = "4215 acMMAE + 4194 unconjugated MMAE concentration-time pairs (Lu 2019 Results section first paragraph)",
    age_range      = "Adults with NHL (Lu 2019 does not tabulate age in the main paper; the four constituent studies enrolled adult patients with NHL or CLL).",
    weight_range   = "5th-95th percentile 48.7-118 kg per Lu 2019 Figure 3 sensitivity-analysis annotation; reference 75 kg used for all weight-based covariate normalizations.",
    sex_female_pct = NA_real_,
    race_ethnicity = "Tested as Asian vs non-Asian indicator only in the final model (RACE_ASIAN); a multiplicative effect on Vc (e_asian_vc = 0.929) was retained based on diagnostic-plot evidence per Lu 2019 Methods step 4 narrative.",
    disease_state  = "Relapsed/refractory or previously untreated B-cell non-Hodgkin lymphoma (NHL): diffuse large B-cell lymphoma (DLBCL) and follicular lymphoma (FL); a small CLL sub-cohort from study DCS4968g is also included.",
    dose_range     = "Pola 0.1-2.4 mg/kg IV every 3 weeks (Q3W) as monotherapy or in combination with rituximab, obinutuzumab, bendamustine, cyclophosphamide, and/or doxorubicin; Q4W in the FL DLBCL+R-bendamustine cohort (study GO29365). Therapeutic dose 1.8 mg/kg Q3W. Dose was administered as a 1.5- to 4-hour infusion.",
    regions        = "Multi-regional (the four constituent studies were global Phase I/Ib/II trials NCT01290549, NCT01691898, NCT02257567, NCT01992653).",
    studies        = "DCS4968g (NCT01290549, Phase I/Ib single-agent and Pola+R run-in), GO27834 / ROMULUS (NCT01691898, Phase Ib/II Pola+R or Pola+G), GO29365 (NCT02257567, Phase Ib/II Pola+B+R or Pola+B+G), GO29044 (NCT01992653, Phase Ib/II Pola+R+CHP or Pola+G+CHP first-line DLBCL). See Lu 2019 Table S1.",
    reference_subject = "75 kg, ALB 35 g/L, TUMSZ 5000 mm^2 SPD, B-cell 1 cell/uL (so max(1, BLBCELL) = 1), male, R/R, non-Asian, normal hepatic function, ECOG >= 1, single-agent. With these covariates COVVC = COVCL_SS = COVKDES = COVCL_TIME = COVMMAE = 1.",
    notes          = "Population characteristics drawn from Lu 2019 main text Results, Figure 3 sensitivity-analysis annotations, and Table S1. The integrated model was developed sequentially (acMMAE base -> integrated acMMAE-MMAE base -> acMMAE covariate -> integrated covariate); the function returns the final integrated covariate model (Tables 1, 2, S3)."
  )

  ini({
    # ----- acMMAE structural parameters (Lu 2019 Table 1, theta1-theta11) -----
    # Reference subject: 75 kg, ALB 35 g/L, TUMSZ 5000 mm^2, B-cell 1 cell/uL,
    # male, R/R, non-Asian, normal hepatic function, ECOG >= 1, single-agent.
    lkdes      <- log(0.0046);  label("Rate constant of CL_TIME exponential decay (kdes, 1/hour)")     # Lu 2019 Table 1, theta1
    lcl_time   <- log(0.00623); label("Initial CL_TIME at time 0 for the reference subject (CL_TIME, L/hour)") # Lu 2019 Table 1, theta2
    lcl     <- log(0.0344);  label("acMMAE nonspecific linear clearance after repeated dosing (CL_SS, L/hour)") # Lu 2019 Table 1, theta3
    lvc        <- log(3.15);    label("acMMAE central volume (Vc, L)")                                 # Lu 2019 Table 1, theta4
    lvp        <- log(3.98);    label("acMMAE peripheral volume (Vp, L)")                              # Lu 2019 Table 1, theta5
    lq         <- log(0.0145);  label("acMMAE intercompartmental clearance (Q, L/hour)")               # Lu 2019 Table 1, theta6
    lvmax      <- log(0.0203);  label("acMMAE Michaelis-Menten maximum elimination rate (Vmax, ng/mL/hour)") # Lu 2019 Table 1, theta7
    lkm_ac     <- log(0.604);   label("acMMAE Michaelis-Menten constant (KM, ng/mL)")                  # Lu 2019 Table 1, theta8
    clss_emax  <- 0.223;        label("Maximum fractional effect of cycle on CL_NS (CLSSEMAX, unitless)") # Lu 2019 Table 1, theta9
    lt50_mo    <- log(3.53);    label("Time of half-maximal cycle effect on CL_NS (T50, months)")      # Lu 2019 Table 1, theta10 (converted to hours inside model() via T50_hr = T50_mo * 24 * 30)
    gamma_ns   <- 2.27;         label("Sigmoidicity of the CL_NS(t) Hill function (gamma, unitless)")  # Lu 2019 Table 1, theta11

    # ----- Unconjugated MMAE structural parameters (Lu 2019 Table 1, theta12-theta21) -----
    # CL_MMAE / V_MMAE / Q_MMAE / V2_MMAE are *apparent* parameters: the
    # absolute fraction of formation of MMAE from acMMAE cannot be estimated,
    # so the systemic CL/V values are scaled by 1 / (true fraction of formation).
    lvc_mmae   <- log(82.2);    label("Unconjugated MMAE apparent central volume (V_MMAE, L)")         # Lu 2019 Table 1, theta12
    lcl_mmae   <- log(1.89);    label("Unconjugated MMAE apparent linear clearance (CL_MMAE, L/hour)") # Lu 2019 Table 1, theta13
    lq_mmae    <- log(36.3);    label("Unconjugated MMAE apparent intercompartmental clearance (Q_MMAE, L/hour)") # Lu 2019 Table 1, theta14
    lvp_mmae   <- log(200);     label("Unconjugated MMAE apparent peripheral volume (V2_MMAE, L)")     # Lu 2019 Table 1, theta15
    lvmax_mmae <- log(0.0307);  label("Unconjugated MMAE Michaelis-Menten maximum elimination rate (Vmax_MMAE, ng/mL/hour)") # Lu 2019 Table 1, theta16
    lkss_mmae  <- log(0.581);   label("Unconjugated MMAE Michaelis-Menten constant (KSS, ng/mL)")      # Lu 2019 Table 1, theta17
    lfrac_clt  <- log(3.70);    label("FRAC_CLT: ratio of acMMAE-MMAE conversion fraction for CL_t pathway relative to CL_NS (unitless)") # Lu 2019 Table 1, theta18
    lfrac_mm   <- log(2.72);    label("FRAC_MM: ratio of acMMAE-MMAE conversion fraction for CL_MM pathway relative to CL_NS (unitless)") # Lu 2019 Table 1, theta19
    lalph_mo   <- log(0.167);   label("Rate constant of FRAC_NS time-decay (alpha, 1/month)")          # Lu 2019 Table 1, theta20 (converted to 1/hour inside model() via alpha_hr = alpha_mo / (24 * 30))
    frac_t     <- 0.139;        label("Initial-time-dependent multiplier of FRAC_NS (FRAC_T, unitless)") # Lu 2019 Table 1, theta21

    # ----- Covariate effects on acMMAE parameters (Lu 2019 Table 2, theta22-theta37) -----
    # WT enters all acMMAE clearance and volume parameters as power effects
    # normalized to 75 kg.
    e_wt_cl      <-  0.73;   label("Power exponent of WT on CL_SS (unitless)")                      # Lu 2019 Table 2, theta22
    e_wt_vc         <-  0.50;   label("Shared power exponent of WT on Vc, Vp, Q (unitless)")           # Lu 2019 Table 2, theta23

    # SEXF coding inverts vs. the source's male-indicator: paper theta24 = 1.20
    # is V1_male / V1_female; with SEXF the multiplicative ratio is V1_female /
    # V1_male = 1 / 1.20 = 0.8333. theta27 = 1.10 -> 1/1.10 = 0.9091. theta39 =
    # 0.911 -> 1/0.911 = 1.0977. (See Errata / Assumptions section in vignette.)
    e_sexf_vc       <-  1 / 1.20;  label("Multiplicative effect of female sex on Vc (ratio female:male, unitless)")     # Lu 2019 Table 2, theta24 inverted (paper reports 1.20 as male:female; SEXF stores 1/1.20)
    e_asian_vc      <-  0.929;     label("Multiplicative effect of Asian race on Vc (unitless)")                        # Lu 2019 Table 2, theta25
    e_line1l_vc     <-  1.20;      label("Multiplicative effect of treatment-naive (first-line) status on Vc (unitless)") # Lu 2019 Table 2, theta26
    e_sexf_cl    <-  1 / 1.10;  label("Multiplicative effect of female sex on CL_SS (ratio female:male, unitless)")  # Lu 2019 Table 2, theta27 inverted
    e_alb_cl     <- -0.247;     label("Power exponent of ALB on CL_SS (unitless; reference 35 g/L)")                 # Lu 2019 Table 2, theta28
    e_combo_rg_cl <- 0.844;     label("Multiplicative effect of anti-CD20 (rituximab or obinutuzumab) combination on CL_SS (unitless)") # Lu 2019 Table 2, theta29
    e_blbcell_cl <-  0.0212;    label("Power exponent of max(1, BLBCELL) on CL_SS (unitless; B-cell count in cells/uL floored at 1)") # Lu 2019 Table 2, theta30
    e_tumsz_cl   <-  0.0521;    label("Linear coefficient of (TUMSZ/5000 - 1) on CL_SS (unitless; effect = 1 + theta * (TUMSZ/5000 - 1), reference 5000 mm^2 SPD)") # Lu 2019 Table 2, theta31
    e_line1l_kdes   <-  3.38;      label("Multiplicative effect of treatment-naive status on kdes (unitless)")          # Lu 2019 Table 2, theta32
    e_combo_rg_kdes <-  0.932;     label("Multiplicative effect of anti-CD20 combination on kdes (unitless)")           # Lu 2019 Table 2, theta33
    e_line1l_cl_time <- 3.53;      label("Multiplicative effect of treatment-naive status on CL_TIME (unitless)")        # Lu 2019 Table 2, theta34
    tmbd50_cl_time  <-  1150;      label("Half-maximal-effect TUMSZ on CL_TIME (Michaelis-Menten-style scaling, mm^2 SPD; effect = TUMSZ / (tmbd50_cl_time + TUMSZ))") # Lu 2019 Table 2, theta35
    bcell_thr_cl_time <- 121;      label("B-cell threshold below which BLBCELL has no effect on CL_TIME (cells/uL; effect = max(1, BLBCELL/threshold)^exponent)") # Lu 2019 Table 2, theta36
    e_blbcell_cl_time <- 0.578;    label("Power exponent of max(1, BLBCELL/threshold) on CL_TIME (unitless)")           # Lu 2019 Table 2, theta37

    # ----- Reference (covariate-free) value of FRAC_NS -----
    # FRAC_NS (the acMMAE -> MMAE conversion fraction) has no estimated typical
    # value in Lu 2019; the reference value is implicitly 1.0 (covariates scale
    # off this reference). The fixed log-parameter lfrac_mmae = log(1) = 0 is
    # carried so the IIV `etalfrac_mmae` has a matching structural parameter
    # per nlmixr2lib convention; setting it to fixed(0) preserves the original
    # behavior (FRAC_0 = cov_mmae * exp(etalfrac_mmae) when lfrac_mmae = 0).
    lfrac_mmae <- fixed(log(1)); label("Reference (covariate-free) FRAC_NS on log scale (unitless; FIXED at 0)")

    # ----- Covariate effects on FRAC_NS (acMMAE -> MMAE conversion fraction; Lu 2019 Table 2, theta38-theta44) -----
    e_wt_frac_mmae        <- -0.467;    label("Power exponent of WT on FRAC_NS (unitless; reference 75 kg)")                  # Lu 2019 Table 2, theta38
    e_sexf_frac_mmae      <-  1 / 0.911; label("Multiplicative effect of female sex on FRAC_NS (ratio female:male, unitless)") # Lu 2019 Table 2, theta39 inverted
    e_line1l_frac_mmae    <-  0.756;    label("Multiplicative effect of treatment-naive status on FRAC_NS (unitless)")        # Lu 2019 Table 2, theta40
    e_combo_rg_frac_mmae  <-  0.709;    label("Multiplicative effect of anti-CD20 combination on FRAC_NS (unitless)")         # Lu 2019 Table 2, theta41
    e_hepimp_frac_mmae    <-  1.19;     label("Multiplicative effect of NCI ODWG hepatic impairment on FRAC_NS (unitless)")   # Lu 2019 Table 2, theta42
    e_ecog_ge1_frac_mmae  <-  1 / 0.905; label("Multiplicative effect of ECOG_GE1 (ECOG >= 1) on FRAC_NS (ratio ECOG>=1:ECOG=0, unitless)") # Lu 2019 Table 2, theta43 inverted
    e_alb_frac_mmae       <- -0.613;    label("Power exponent of ALB on FRAC_NS (unitless; reference 35 g/L)")                # Lu 2019 Table 2, theta44

    # ----- Inter-individual variability (Lu 2019 Table S3, omega^2 values) -----
    # Stored as variances on the log scale exactly as the paper reports them
    # (the paper labels these omega^2; the per-row %CV column matches sqrt(exp(omega^2) - 1)).
    etalcl_time ~ 1.89                       # Lu 2019 Table S3, Omega11 (CV 138%)
    etalcl  ~ 0.0376                      # Lu 2019 Table S3, Omega22 (CV 19.5%)
    etalvc     ~ 0.0151                      # Lu 2019 Table S3, Omega33 (CV 12.3%)
    etalvp     ~ 0.107                       # Lu 2019 Table S3, Omega44 (CV 32.7%)
    etalq      ~ 0.0538                      # Lu 2019 Table S3, Omega55 (CV 23.2%)
    etalvmax   ~ 0.462                       # Lu 2019 Table S3, Omega66 (CV 67.9%)
    etalfrac_mmae ~ 0.0972                   # Lu 2019 Table S3, Omega77 -- IIV on FRAC_0 = COVMMAE * exp(eta) (CV 31.2%)
    etalcl_mmae ~ 0.115                      # Lu 2019 Table S3, Omega88 (CV 33.9%)
    etalvp_mmae ~ 0.0422                     # Lu 2019 Table S3, Omega99 (CV 20.5%)

    # ----- Residual error -----
    # Lu 2019 NONMEM uses an eta-on-epsilon parameterization: Y = TY * (1 +
    # EPS * exp(ETA10 or ETA11)). The two residual epsilons are correlated via
    # the omega(10,10), omega(10,11), omega(11,11) block of Table S3. nlmixr2
    # does not natively support eta-on-epsilon; we collapse to a fixed
    # proportional residual error per analyte using the Sigma point estimates
    # (sqrt(Sigma)) and document this deviation in the validation vignette.
    propSd       <- sqrt(0.0254); label("Proportional residual error on acMMAE Cc (fraction)")             # Lu 2019 Table S3, Sigma11 = 0.0254
    propSd_mmae    <- sqrt(0.0726); label("Proportional residual error on unconjugated MMAE Cc_mmae (fraction)") # Lu 2019 Table S3, Sigma22 = 0.0726
  })

  model({
    # ===== 1. Derived covariate terms (Lu 2019 supplement, NONMEM $PK / Notations block) =====
    # B-cell flooring: BCEL is the BLBCELL/threshold ratio floored at 1 (used
    # in the CL_TIME covariate); BCEL1 is BLBCELL itself floored at 1 cell/uL
    # (used in the CL_SS covariate). max(1, x) preserves the NONMEM
    # IF(BBCC.GT.thresh) BCEL=BBCC/thresh; ELSE BCEL=1 logic.
    bcel_cl_time <- max(1, BLBCELL / bcell_thr_cl_time)  # Lu 2019 supplement: BCEL  for CL_TIME
    bcel_cl   <- max(1, BLBCELL)                      # Lu 2019 supplement: BCEL1 for CL_SS

    # Aggregate covariate multipliers (Lu 2019 supplement, COVV1, COVCLINF,
    # COVKDES, COVCLT, COVMMAE expressions), with reference categories at 75 kg
    # / 35 g/L / 5000 mm^2 SPD / 1 cell/uL / male / R/R / non-Asian / ECOG>=1 /
    # normal hepatic function / single agent.
    cov_vc <- (WT / 75)^e_wt_vc *
              e_sexf_vc^SEXF *
              e_asian_vc^RACE_ASIAN *
              e_line1l_vc^LINE_1L

    cov_cl <- (WT / 75)^e_wt_cl *
                 e_sexf_cl^SEXF *
                 (ALB / 35)^e_alb_cl *
                 e_combo_rg_cl^COMBO_RG *
                 bcel_cl^e_blbcell_cl *
                 (1 + e_tumsz_cl * (TUMSZ / 5000 - 1))

    cov_kdes <- e_line1l_kdes^LINE_1L *
                e_combo_rg_kdes^COMBO_RG

    cov_cl_time <- e_line1l_cl_time^LINE_1L *
                   (TUMSZ / (tmbd50_cl_time + TUMSZ)) *
                   bcel_cl_time^e_blbcell_cl_time

    cov_mmae <- (WT / 75)^e_wt_frac_mmae *
                e_sexf_frac_mmae^SEXF *
                e_line1l_frac_mmae^LINE_1L *
                e_combo_rg_frac_mmae^COMBO_RG *
                e_hepimp_frac_mmae^HEPIMP *
                e_ecog_ge1_frac_mmae^ECOG_GE1 *
                (ALB / 35)^e_alb_frac_mmae

    # ===== 2. Individual PK parameters (Lu 2019 supplement, Random Effects Model block) =====
    # acMMAE side
    kdes        <- exp(lkdes) * cov_kdes                      # 1/hour
    cl_time_init <- exp(lcl_time + etalcl_time) * cov_cl_time # CL_TIME initial value at t = 0 (L/hour)
    cl       <- exp(lcl + etalcl) * cov_cl        # L/hour
    vc          <- exp(lvc + etalvc) * cov_vc                 # L
    vp          <- exp(lvp + etalvp) * (WT / 75)^e_wt_vc      # L
    q           <- exp(lq + etalq) * (WT / 75)^e_wt_vc        # L/hour
    vmax        <- exp(lvmax + etalvmax)                      # ng/mL/hour
    km_ac       <- exp(lkm_ac)                                # ng/mL
    t50         <- exp(lt50_mo) * 24 * 30                     # T50 in hours (paper reports months)
    t50gam      <- t50^gamma_ns

    # Unconjugated MMAE side (apparent parameters)
    vc_mmae   <- exp(lvc_mmae)                            # L
    cl_mmae   <- exp(lcl_mmae + etalcl_mmae)              # L/hour
    q_mmae    <- exp(lq_mmae)                             # L/hour
    vp_mmae   <- exp(lvp_mmae + etalvp_mmae)              # L
    vmax_mmae <- exp(lvmax_mmae)                          # ng/mL/hour
    kss_mmae  <- exp(lkss_mmae)                           # ng/mL
    frac_clt  <- exp(lfrac_clt)                           # unitless
    frac_mm   <- exp(lfrac_mm)                            # unitless
    alph      <- exp(lalph_mo) / (24 * 30)                # alpha in 1/hour (paper reports 1/month)
    frac_0    <- exp(lfrac_mmae + etalfrac_mmae) * cov_mmae  # FRAC_0 = COVMMAE * exp(eta7); lfrac_mmae fixed at 0 keeps original behavior

    # ===== 3. Time-dependent quantities (Lu 2019 supplement $DES block) =====
    # Hill function for CL_NS(t) (Lu 2019 Eq. 1). At t = 0, CL_NS = CL_SS *
    # (1 + CLSSEMAX); as t -> infinity, CL_NS -> CL_SS.
    tgam <- time^gamma_ns
    cl_ns <- cl * (1 + clss_emax * t50gam / (t50gam + tgam))

    # Exponential decay of CL_TIME (Lu 2019 supplement Notations: CLT = CLT0 * exp(-kdes*t)).
    cl_t <- cl_time_init * exp(-kdes * time)

    # Time-dependent FRAC_NS (Lu 2019 supplement: FRAC = FRAC_0 * (1 + FRAC_T * exp(-alpha*t))).
    frac_ns <- frac_0 * (1 + frac_t * exp(-alph * time))

    # ===== 4. Micro-constants =====
    k12 <- q / vc
    k21 <- q / vp
    k34 <- q_mmae / vc_mmae
    k43 <- q_mmae / vp_mmae
    k30 <- cl_mmae / vc_mmae
    k10 <- (cl_t + cl_ns) / vc                            # Linear (non-MM) elimination of acMMAE

    # KINPUT = FRAC_NS * (CL_NS + FRAC_CLT * CL_TIME + FRAC_MM * CL_MM) / Vc
    # where CL_MM = Vmax * Vc / (KM + A1/Vc).
    # The product KINPUT * A1 is the input rate (amount/hour) into central MMAE
    # contributed by the three acMMAE elimination pathways (Lu 2019 Eq. 2).
    kinput <- frac_ns * (
                cl_ns +
                frac_clt * cl_t +
                frac_mm * vmax * vc / (km_ac + central / vc)
              ) / vc

    # ===== 5. ODE system (Lu 2019 supplement $DES; signs verified against NONMEM DADT lines) =====
    # NOTE: the published "Equations" section of the supplement contains a sign
    # typo on the dA2/dt line ("dA2/dt = -K12 A1 - K21 A2" should read
    # "dA2/dt = +K12 A1 - K21 A2"). The NONMEM control stream's DADT(2) line
    # (DADT(2) = -K21*A(2) + K12*A(1)) gives the correct mass-conservation
    # form, which we adopt here. See the Errata section of the validation
    # vignette for details.
    d/dt(central)    <- -k10 * central -
                         k12 * central +
                         k21 * peripheral1 -
                         vmax * central / (km_ac + central / vc)
    d/dt(peripheral1) <-  k12 * central -
                          k21 * peripheral1

    d/dt(central_mmae)    <-  kinput * central -
                              k30 * central_mmae -
                              k34 * central_mmae +
                              k43 * peripheral1_mmae -
                              vmax_mmae * central_mmae / (kss_mmae + central_mmae / vc_mmae)
    d/dt(peripheral1_mmae) <- k34 * central_mmae -
                              k43 * peripheral1_mmae

    # ===== 6. Observations and residual error =====
    # Concentration: amount in ug, volume in L -> ug/L = ng/mL.
    Cc      <- central / vc
    Cc_mmae <- central_mmae / vc_mmae

    Cc      ~ prop(propSd)
    Cc_mmae ~ prop(propSd_mmae)
  })
}
