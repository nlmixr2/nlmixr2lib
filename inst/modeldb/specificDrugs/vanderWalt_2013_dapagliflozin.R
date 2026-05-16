vanderWalt_2013_dapagliflozin <- function() {
  description <- paste(
    "Semi-mechanistic joint parent-metabolite population PK model for",
    "dapagliflozin and its inactive UGT1A9 glucuronide metabolite",
    "dapagliflozin 3-O-glucuronide (D3OG, identified as M15 in",
    "chromatography) in healthy adults, T2DM subjects with normal or",
    "impaired renal function, and patients with hepatic impairment",
    "(van der Walt 2013). Parent: 2-compartment disposition with",
    "first-order absorption fed by a Savic 2007 transit-compartment",
    "chain (continuous N estimated alongside MTT) and a logit",
    "bioavailability anchor; three parallel parent elimination",
    "pathways are estimated separately as renal excretion of unchanged",
    "dapagliflozin (CLP_renal, proportional to baseline creatinine",
    "clearance), metabolic formation of D3OG (CLP_M15), and metabolic",
    "clearance to unmeasured metabolites (CLP_other, allometrically",
    "scaled like CLP_M15). Metabolite: 1-compartment with renal",
    "elimination CLM proportional to creatinine clearance. Plasma",
    "observations only are emulated here -- the source paper also",
    "fitted urine dapagliflozin and D3OG concentrations simultaneously",
    "with a replicate residual-error structure; see the validation",
    "vignette for the urine and replicate-residual deviations.",
    "Covariates: creatinine clearance (CRCL; IBW-corrected, mL/min)",
    "on CLP_M15, CLP_renal, and CLM; AGE on CLP_other;",
    "Child-Pugh Class C (HEPIMP_SEV) on CLP_M15 and V2M;",
    "Child-Pugh Class B or C (HEPIMP_MODSEV) on V3P and CLM; female",
    "sex (SEXF) on total CLP and on CLM; allometric WT scaling on",
    "CLP_M15, CLP_other, V2P, V3P, V2M.",
    sep = " "
  )
  reference <- paste(
    "van der Walt J-S, Hong Y, Zhang L, Pfister M, Boulton DW, Karlsson MO.",
    "A nonlinear mixed effects pharmacokinetic model for dapagliflozin",
    "and dapagliflozin 3O-glucuronide in renal or hepatic impairment.",
    "CPT Pharmacometrics Syst Pharmacol. 2013;2(5):e42.",
    "doi:10.1038/psp.2013.20.",
    sep = " "
  )
  vignette <- "vanderWalt_2013_dapagliflozin"
  units <- list(
    time          = "hour",
    dosing        = "mg",
    concentration = "ng/mL"
  )

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed baseline body weight. Used for a priori allometric",
        "scaling on metabolic clearances CLP_M15 and CLP_other with",
        "exponent 0.75 and on volumes V2P, V3P, and V2M with exponent",
        "1.0, all referenced to 70 kg (van der Walt 2013 Methods text",
        "'Metabolic clearances were scaled by the individual baseline",
        "body weight using an allometric model, (WT/70)' and Table 2",
        "footnote a 'Allometric scaling (BBWT_i/70)^(3/4)'). CLP_renal",
        "and CLM are not scaled by body weight (Table 2: 'x' on the",
        "WT column for CLP_renal; CLM is driven only by CL_cr_IBW)."
      ),
      source_name        = "BBWT"
    ),
    AGE = list(
      description        = "Baseline age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed baseline age. Linear additive effect on CLP_other",
        "of (-0.0204) per year of difference from the reference age",
        "53.5 years (van der Walt 2013 Table 2 reference value;",
        "paper text 'CLP_other decreased by 2% for every year of",
        "increasing age above 54 years'). The covariate effect is",
        "applied as (1 + e_age_cl_nonren * (AGE - 53.5))."
      ),
      source_name        = "AGE"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male) per the canonical SEXF column. van der Walt 2013 Table 2 also uses male as the reference category and reports the female multiplicative shift.",
      notes              = paste(
        "Time-fixed. Multiplicative fractional effects on total",
        "dapagliflozin clearance CLP (applied uniformly to each",
        "component CLP_renal / CLP_M15 / CLP_other so that the sum",
        "shifts by the paper-reported factor) and on D3OG renal",
        "clearance CLM. In females (SEXF = 1) CLP is 16.7% lower",
        "(e_female_cl = -0.167) and CLM is 19.6% lower",
        "(e_female_cl_d3og = -0.196); van der Walt 2013 Final model",
        "(Table 1) and Results 'In females, CLP and CLM were 16.7 and",
        "19.6% lower, respectively'. Reference is male (SEXF = 0)."
      ),
      source_name        = "SEX"
    ),
    CRCL = list(
      description        = "Baseline creatinine clearance, ideal-body-weight-corrected (Cockcroft-Gault with IBW)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed baseline value. Computed by van der Walt 2013 using",
        "the Cockcroft-Gault equation with ideal body weight (IBW) in",
        "place of total body weight (BBWT); reported as 'CL_cr,IBW' in",
        "the source paper and supplied as the NONMEM column for the",
        "creatinine-clearance covariate. Reported units are mL/min --",
        "NOT BSA-normalized (i.e., values are not divided by BSA/1.73).",
        "Used with three different functional forms: (i) linear in CRCL",
        "with intercept zero on CLP_renal (cl_renal = th_cl_renal *",
        "CRCL, where th_cl_renal is in (L/h)/(mL/min)); (ii) linear in",
        "CRCL with intercept zero on CLM (cl_d3og = th_cl_d3og * CRCL,",
        "similarly); (iii) proportional-change-from-the-mean (centered",
        "at 80.14 mL/min, the pooled-cohort mean) on the dapagliflozin",
        "to D3OG formation clearance CLP_M15, applied as (1 + e_crcl *",
        "(CRCL - 80.14)). Population range across the three contributing",
        "studies: 13 to 143 mL/min."
      ),
      source_name        = "CL_cr,IBW"
    ),
    HEPIMP_SEV = list(
      description        = "Severe hepatic impairment indicator (Child-Pugh Class C, score 10-15)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (any non-severe Child-Pugh category: normal, mild = Child-Pugh A, or moderate = Child-Pugh B). The HEPIMP_SEV = 1 cohort in the source study is the six subjects with Child-Pugh Class C in the hepatic-impairment substudy.",
      notes              = paste(
        "Time-fixed binary indicator of severe hepatic impairment per",
        "the Child-Pugh classification (Class C). NOT the NCI ODWG",
        "default classification carried by the canonical HEPIMP_SEV",
        "column (see inst/references/covariate-columns.md HEPIMP_SEV",
        "entry, source aliases). Multiplicative fractional effects:",
        "CLP_M15 decreased by 42.2% in HEPIMP_SEV = 1 cohort",
        "(e_hepsev_cl_form_d3og = -0.422) and V2M increased by 134%",
        "(e_hepsev_vc_d3og = +1.33); van der Walt 2013 Final model",
        "(Table 1, Covariate fixed effects) and Results 'With severe",
        "HI (Child-Pugh Class C), CLP_M15 decreased by 41% and V2M",
        "increased by 134%'. Note the per-model semantics: HEPIMP_SEV",
        "is paired with HEPIMP_MODSEV (Child-Pugh Class B or C) so",
        "that HEPIMP_SEV = 1 ALWAYS implies HEPIMP_MODSEV = 1 (a",
        "Child-Pugh C subject is also Child-Pugh B-or-C). Reference",
        "category for both indicators is 'not Child-Pugh C and not",
        "Child-Pugh B-or-C' = normal or mild hepatic function."
      ),
      source_name        = "Child-Pugh Class C"
    ),
    HEPIMP_MODSEV = list(
      description        = "Composite moderate-or-severe hepatic impairment indicator (Child-Pugh Class B or C, score >= 7)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (normal hepatic function or mild Child-Pugh Class A). The HEPIMP_MODSEV = 1 cohort in the source study pools 12 subjects (6 Child-Pugh B + 6 Child-Pugh C) of the 24 subjects in the hepatic-impairment substudy.",
      notes              = paste(
        "Time-fixed binary indicator of moderate-or-severe hepatic",
        "impairment per the Child-Pugh classification (Class B or C",
        "pooled). NOT the NCI ODWG default classification carried by",
        "the canonical HEPIMP_MODSEV column (see",
        "inst/references/covariate-columns.md HEPIMP_MODSEV entry,",
        "source aliases). Multiplicative fractional effects: V3P",
        "decreased by 60.0% in HEPIMP_MODSEV = 1 cohort",
        "(e_hepmodsev_vp = -0.600) and CLM decreased by 29.3%",
        "(e_hepmodsev_cl_d3og = -0.293); van der Walt 2013 Final model",
        "(Table 1, Covariate fixed effects) and Results 'Moderate or",
        "severe HI (Child-Pugh Class B or C) decreased CLM and the",
        "peripheral volume of distribution of dapagliflozin (V3P) by",
        "29 and 60%, respectively'. Paired with HEPIMP_SEV: a",
        "Child-Pugh C subject carries HEPIMP_SEV = 1 AND",
        "HEPIMP_MODSEV = 1; a Child-Pugh B subject carries HEPIMP_SEV",
        "= 0 AND HEPIMP_MODSEV = 1. The two indicators are nested,",
        "not mutually exclusive."
      ),
      source_name        = "Child-Pugh Class B or C"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 227L,
    n_studies      = 3L,
    age_range      = "25-92 years (per-study medians 63, 43, 67)",
    age_median     = "67 years (largest study, MB102029)",
    weight_range   = "51.8-148.3 kg",
    weight_median  = "approximately 86 kg (pooled across the three studies; per-study medians 81.2, 86.3, 92.0 kg)",
    sex_female_pct = 32.6,
    race_ethnicity = NULL,
    disease_state  = paste(
      "Pooled cohort of (i) healthy adult subjects, (ii) adults with",
      "type 2 diabetes mellitus (T2DM) with normal, mild, moderate, or",
      "severe renal impairment, and (iii) adults with hepatic",
      "impairment matched to healthy controls. Study MB102007:",
      "NCT00554450 renal-impairment PK/PD study, 40 subjects (8",
      "healthy + 32 T2DM at varying renal function), single 50-mg dose",
      "then 20-mg q.d. for 7 days. Study MB102027: 5-day single-dose",
      "10-mg hepatic-impairment study, 24 subjects (18 HI + 6 healthy",
      "matched controls), 12 with mild/moderate impairment and 12",
      "with severe impairment per Child-Pugh classification. Study",
      "MB102029: NCT00663260 52-week phase 3 trial in 163 T2DM",
      "subjects with moderate renal impairment."
    ),
    renal_function = "Pooled cohort spans CRCL 13-143 mL/min (Table 2 baseline-covariate summary).",
    hepatic_function = "MB102027 substudy: 12 normal-or-Pugh-A, 6 Pugh B, 6 Pugh C; MB102007 and MB102029 subjects were normal or mild only.",
    dose_range     = paste(
      "Oral dapagliflozin: 50-mg single dose followed by 20-mg q.d.",
      "for 7 days (renal-impairment PK/PD study MB102007); 10-mg",
      "single dose (hepatic-impairment study MB102027); 10-mg daily",
      "(phase 3 trial MB102029)."
    ),
    regions        = "Not specified by region in the publication; multi-centre Bristol-Myers Squibb clinical program.",
    notes          = paste(
      "Demographics from van der Walt 2013 Table 2 baseline-covariate",
      "summary across the three contributing studies. NONMEM 7.1.2",
      "with FOCE and eta-eps interaction; standard errors via Monte",
      "Carlo importance sampling EM assisted by MAP estimation",
      "(IMPMAP) with 10,000 samples per individual. Note: nlmixr2lib",
      "extraction emulates plasma observations only (Cc for",
      "dapagliflozin, Cc_d3og for D3OG); urine excretion modeling",
      "and the cross-output replicate residual-error structure",
      "reported in the paper are out of scope for this extraction --",
      "see vignette Errata."
    )
  )

  ini({
    # ============================================================
    # Parent dapagliflozin structural parameters -- van der Walt
    # 2013 Table 1, "Final model" typical-value column. Reference
    # subject for typical-value parameters: 70-kg male (SEXF = 0),
    # age 53.5 years, CRCL = 80.14 mL/min, normal hepatic function
    # (HEPIMP_SEV = HEPIMP_MODSEV = 0).
    # ============================================================

    # Renal clearance of unchanged parent -- linear in CRCL with
    # intercept zero (so cl_renal = th_cl_renal * CRCL). Reported in
    # (L/h)/(mL/min) so that the same CRCL column in mL/min drives
    # the L/h clearance directly.
    lth_cl_renal <- log(0.00310)
    label("Parent renal-clearance coefficient (L/h per mL/min CRCL)")  # Table 1 Final model: CLP_renal = 0.00310 (L/h)/(mL/min) RSE 9.9%

    # Parent metabolic clearance to D3OG (formation), at the 70-kg
    # reference body weight and at the 80.14 mL/min reference CRCL.
    lcl_form_d3og <- log(7.54)
    label("Parent -> D3OG formation clearance CLP_M15 at WT = 70 kg, CRCL = 80.14 mL/min, normal hepatic function (L/h)")  # Table 1 Final model: CLP_M15 = 7.54 L/h RSE 12.0%

    # Parent metabolic clearance to unmeasured metabolites, at the
    # 70-kg reference body weight and at the 53.5-year reference age.
    lcl_nonren <- log(5.35)
    label("Parent metabolic clearance to unmeasured metabolites CLP_other at WT = 70 kg, age = 53.5 years (L/h)")  # Table 1 Final model: CLP_other = 5.35 L/h RSE 12.1%

    # Parent central and peripheral disposition.
    lvc <- log(39.0)
    label("Parent central volume of distribution V2P at WT = 70 kg (L)")  # Table 1 Final model: V2P = 39.0 L RSE 4.4%
    lq <- log(7.07)
    label("Parent inter-compartmental clearance QP (L/h)")  # Table 1 Final model: QP = 7.07 L/h RSE 7.1%
    lvp <- log(71.5)
    label("Parent peripheral volume of distribution V3P at WT = 70 kg, normal-or-mild hepatic function (L)")  # Table 1 Final model: V3P = 71.5 L RSE 7.0%

    # Absorption: Savic 2007 transit-compartment chain. MTT is the
    # mean transit time and NN is the (continuous) number of transit
    # compartments preceding first-order absorption into central.
    lmtt <- log(0.475)
    label("Mean transit time MTT through the absorption chain (h)")  # Table 1 Final model: MTT = 0.475 h RSE 7.2%
    lnn <- log(5.45)
    label("Number of transit compartments NN (continuous, dimensionless)")  # Table 1 Final model: NN = 5.45 RSE 16.5%

    # Bioavailability: logit-transformed so the back-transformed F
    # stays in (0, 1). Typical value BIO = 0.858 (final model). The
    # logit-scale typical value is log(0.858 / (1 - 0.858)) = 1.799.
    logitfdepot <- log(0.858 / (1 - 0.858))
    label("Logit-transformed typical bioavailability F (logit-scale; back-transformed F = 0.858)")  # Table 1 Final model: BIO = 0.858 (logit Eq. 6, footnote g)

    # ============================================================
    # D3OG (metabolite) structural parameters -- van der Walt
    # 2013 Table 1, Final model column. The D3OG state carries
    # dapagliflozin mass-equivalents internally (1:1 molar transfer
    # from CLP_M15 flux); V2M and CLM absorb the MW ratio
    # implicitly per the paper's NONMEM convention, so concentrations
    # output in ng/mL match the paper's plasma D3OG observations
    # without an additional MW conversion factor in the observation
    # equation.
    # ============================================================
    lth_cl_d3og <- log(0.0799)
    label("D3OG renal-clearance coefficient (L/h per mL/min CRCL)")  # Table 1 Final model: CLM = 0.0799 (L/h)/(mL/min) RSE 7.6%
    lvc_d3og <- log(2.26)
    label("D3OG central volume of distribution V2M at WT = 70 kg, normal hepatic function (L)")  # Table 1 Final model: V2M = 2.26 L RSE 5.6%

    # ============================================================
    # Covariate fixed effects -- van der Walt 2013 Table 1, Final
    # model column. Continuous-covariate effects are proportional-
    # change-from-the-mean (Table 2 footnote e); categorical-
    # covariate effects are proportional-change-from-reference
    # (footnote f).
    # ============================================================
    e_crcl_cl_form_d3og <- 0.00502
    label("Proportional-change coefficient: CRCL on CLP_M15, centered at 80.14 mL/min (per mL/min)")  # Table 1 Final model: CLP_M15-Baseline CL_cr_IBW = 0.00502 RSE 42.6%
    e_hepsev_cl_form_d3og <- -0.422
    label("Multiplicative fractional shift on CLP_M15 for HEPIMP_SEV = 1 (unitless)")  # Table 1 Final model: CLP_M15-Child-Pugh Class C = -0.422 RSE 38.6%
    e_age_cl_nonren <- -0.0204
    label("Linear-additive coefficient: AGE on CLP_other, centered at 53.5 years (per year)")  # Table 1 Final model: CLP_other-age = -0.0204 RSE 22.1%
    e_hepsev_vc_d3og <- 1.33
    label("Multiplicative fractional shift on V2M for HEPIMP_SEV = 1 (unitless)")  # Table 1 Final model: V2M-Child-Pugh Class C = 1.33 RSE 17.3%
    e_hepmodsev_vp <- -0.600
    label("Multiplicative fractional shift on V3P for HEPIMP_MODSEV = 1 (unitless)")  # Table 1 Final model: V3P-Child-Pugh Class B,C = -0.600 RSE 12.6%
    e_hepmodsev_cl_d3og <- -0.293
    label("Multiplicative fractional shift on CLM for HEPIMP_MODSEV = 1 (unitless)")  # Table 1 Final model: CLM-Child-Pugh Class B,C = -0.293 RSE 23.3%
    e_female_cl <- -0.167
    label("Multiplicative fractional shift on total parent clearance CLP for SEXF = 1 (unitless)")  # Table 1 Final model: CLP-gender = -0.167 RSE 64.1%
    e_female_cl_d3og <- -0.196
    label("Multiplicative fractional shift on D3OG renal clearance CLM for SEXF = 1 (unitless)")  # Table 1 Final model: CLM-gender = -0.196 RSE 38.9%

    # Allometric exponents on parent metabolic clearances and on all
    # three volumes. Table 2 footnote a fixes the values per
    # convention: 0.75 on CL-like, 1.0 on V-like. The paper does not
    # report uncertainty on either exponent (a priori scaling).
    e_wt_cl_form_d3og <- fixed(0.75)
    label("Allometric exponent on CLP_M15 (fixed at 3/4 per Table 2 footnote a)")  # Table 2 footnote a: a priori (BBWT/70)^(3/4)
    e_wt_cl_nonren <- fixed(0.75)
    label("Allometric exponent on CLP_other (fixed at 3/4 per Table 2 footnote a)")  # Table 2 footnote a: a priori (BBWT/70)^(3/4)
    e_wt_vc <- fixed(1.0)
    label("Allometric exponent on V2P (fixed at 1 per Table 2 a priori scaling)")  # Table 2: a priori (BBWT/70)^1 on V2P
    e_wt_vp <- fixed(1.0)
    label("Allometric exponent on V3P (fixed at 1 per Table 2 a priori scaling)")  # Table 2: a priori (BBWT/70)^1 on V3P
    e_wt_vc_d3og <- fixed(1.0)
    label("Allometric exponent on V2M (fixed at 1 per Table 2 a priori scaling)")  # Table 2: a priori (BBWT/70)^1 on V2M

    # ============================================================
    # IIV -- van der Walt 2013 Table 1, Final model column. Reported
    # as %CV (= sqrt(exp(omega^2) - 1) for log-normal); the internal
    # variance is omega^2 = log(1 + CV^2).
    # ============================================================
    # CLP_M15 and V2M carry a correlation r = 0.562 (Table 1 row
    # "Correlation CLP_M15 : V2M"). Block lower-triangle order in
    # nlmixr2 c(...) is (var_1, cov, var_2).
    #   var_form_d3og = log(1 + 0.368^2) = 0.1270
    #   var_vc_d3og   = log(1 + 0.447^2) = 0.1822
    #   cov           = 0.562 * sqrt(0.1270 * 0.1822) = 0.0855
    etalcl_form_d3og + etalvc_d3og ~ c(0.1270, 0.0855, 0.1822)  # Final model: IIV.CLP_M15 = 36.8% CV, IIV.V2M = 44.7% CV, correlation r = 0.562

    etalcl_nonren ~ 0.0935    # Final model: IIV.CLP_other = 31.3% CV -> omega^2 = log(1 + 0.313^2) = 0.0935
    etalth_cl_renal ~ 0.3243  # Final model: IIV.CLP_renal = 61.9% CV -> omega^2 = log(1 + 0.619^2) = 0.3243
    etalvc ~ 0.0542           # Final model: IIV.V2P = 23.6% CV -> omega^2 = log(1 + 0.236^2) = 0.0542
    etalq ~ 0.1423            # Final model: IIV.QP = 39.1% CV -> omega^2 = log(1 + 0.391^2) = 0.1423
    etalvp ~ 0.1942           # Final model: IIV.V3P = 46.3% CV -> omega^2 = log(1 + 0.463^2) = 0.1942
    etalmtt ~ 0.3100          # Final model: IIV.MTT = 60.3% CV -> omega^2 = log(1 + 0.603^2) = 0.3100
    etalnn ~ 1.520            # Final model: IIV.NN = 189% CV -> omega^2 = log(1 + 1.89^2) = 1.520
    etalth_cl_d3og ~ 0.0602   # Final model: IIV.CLM = 24.9% CV -> omega^2 = log(1 + 0.249^2) = 0.0602

    # IIV on bioavailability is reported by the paper as 11.1% (Table
    # 1 footnote g). The omega is on the logit scale; the back-
    # transformed CV is computed by the delta-method approximation
    # CV ~= sqrt(omega) * (1 - F_TV), so
    # omega ~= (0.111 / (1 - 0.858))^2 = (0.111 / 0.142)^2 = 0.611.
    etalogitfdepot ~ 0.611    # Final model: IIV.BIO = 11.1% back-transformed; raw omega on logit scale = 0.611 per footnote g

    # ============================================================
    # Residual error -- van der Walt 2013 Table 1, Final model
    # column. Plasma-only observations are extracted here (urine
    # observations are out of scope; see vignette Errata). The paper
    # also estimated a "replicate" residual cross-correlation between
    # parent and metabolite observations at the same sampling time;
    # this cross-correlation is NOT representable in the standard
    # nlmixr2 add()/prop() syntax and is omitted in this extraction
    # (the per-output proportional and additive components are
    # retained at the paper-reported magnitudes).
    # ============================================================
    propSd <- 0.207
    label("Parent plasma proportional residual SD (fraction)")  # Table 1 Final model: Prop, dapa_plasma = 0.207 RSE 6.4%
    addSd <- 0.465
    label("Parent plasma additive residual SD (ng/mL)")  # Table 1 Final model: Add, dapa_plasma = 0.465 ng/mL RSE 29.9%
    propSd_d3og <- 0.195
    label("D3OG plasma proportional residual SD (fraction)")  # Table 1 Final model: Prop, M15_plasma = 0.195 RSE 6.4%
    addSd_d3og <- 0.585
    label("D3OG plasma additive residual SD (ng/mL)")  # Table 1 Final model: Add, M15_plasma = 0.585 ng/mL RSE 51.3%
  })

  model({
    # ------------------------------------------------------------
    # Body-weight ratio (reference 70 kg).
    # ------------------------------------------------------------
    wt_ratio <- WT / 70

    # ------------------------------------------------------------
    # Covariate multipliers. Continuous-covariate effects are
    # proportional-change-from-the-mean (Table 2 footnote e):
    #   CRCL on CLP_M15 centered at 80.14 mL/min
    #   AGE on CLP_other centered at 53.5 years
    # Categorical-covariate effects are proportional-change-from-
    # reference (Table 2 footnote f).
    # ------------------------------------------------------------
    crcl_eff_cl_form_d3og <- 1 + e_crcl_cl_form_d3og * (CRCL - 80.14)
    hepsev_eff_cl_form_d3og <- 1 + e_hepsev_cl_form_d3og * HEPIMP_SEV
    age_eff_cl_nonren <- 1 + e_age_cl_nonren * (AGE - 53.5)
    hepsev_eff_vc_d3og <- 1 + e_hepsev_vc_d3og * HEPIMP_SEV
    hepmodsev_eff_vp <- 1 + e_hepmodsev_vp * HEPIMP_MODSEV
    hepmodsev_eff_cl_d3og <- 1 + e_hepmodsev_cl_d3og * HEPIMP_MODSEV
    female_eff_cl <- 1 + e_female_cl * SEXF
    female_eff_cl_d3og <- 1 + e_female_cl_d3og * SEXF

    # ------------------------------------------------------------
    # Individual structural parameters. Each parent-clearance
    # component carries its own IIV and the uniform female-sex
    # multiplicative shift that the paper applies to "total CLP".
    # ------------------------------------------------------------
    cl_renal <- exp(lth_cl_renal + etalth_cl_renal) * CRCL * female_eff_cl
    cl_form_d3og <- exp(lcl_form_d3og + etalcl_form_d3og) *
      wt_ratio^e_wt_cl_form_d3og *
      crcl_eff_cl_form_d3og *
      hepsev_eff_cl_form_d3og *
      female_eff_cl
    cl_nonren <- exp(lcl_nonren + etalcl_nonren) *
      wt_ratio^e_wt_cl_nonren *
      age_eff_cl_nonren *
      female_eff_cl

    vc <- exp(lvc + etalvc) * wt_ratio^e_wt_vc
    vp <- exp(lvp + etalvp) * wt_ratio^e_wt_vp * hepmodsev_eff_vp
    q <- exp(lq + etalq)

    cl_d3og <- exp(lth_cl_d3og + etalth_cl_d3og) * CRCL *
      hepmodsev_eff_cl_d3og *
      female_eff_cl_d3og
    vc_d3og <- exp(lvc_d3og + etalvc_d3og) *
      wt_ratio^e_wt_vc_d3og *
      hepsev_eff_vc_d3og

    mtt <- exp(lmtt + etalmtt)
    nn <- exp(lnn + etalnn)
    # Put the mu-referenced sum on its own simple line first so the
    # rxode2 parser recognises the IIV on logitfdepot before the
    # logit back-transform; nlmixr2 otherwise warns that the eta is
    # non-mu-referenced because the eta sits inside exp(-...).
    logit_fdepot_i <- logitfdepot + etalogitfdepot
    fdepot <- 1 / (1 + exp(-logit_fdepot_i))

    # ka collapses the rxode2 transit() output into central without
    # introducing an extra first-order absorption phase that the
    # source paper does not estimate. The Savic 2007 closed form
    # (van der Walt 2013 Eq. 5, k_TRANSIT = (N + 1)/MTT) already
    # encodes the entire absorption process via the gamma-PDF input
    # rate. Following the Wilkins 2008 rifampicin / Svensson 2018
    # rifampicin pattern, transit() feeds depot and depot absorbs
    # into central at rate ka; with ka >> 1/MTT the depot's
    # exponential tail is negligible and the central-input rate
    # tracks the transit() gamma-PDF directly. We set ka = 60 / h
    # (t1/2 in depot ~= 0.012 h = 0.7 min, an order of magnitude
    # faster than MTT = 0.475 h), well within numerical resolution.
    ka <- 60

    # ------------------------------------------------------------
    # ODEs. Parent: 2-compartment with transit-fed first-order
    # absorption. Metabolite D3OG: 1-compartment fed by the parent's
    # CLP_M15 flux (cl_form_d3og * central / vc). Total parent
    # elimination = (cl_renal + cl_form_d3og + cl_nonren) * central
    # / vc; only the cl_form_d3og component delivers mass into the
    # D3OG compartment (cl_renal and cl_nonren leave the system).
    # ------------------------------------------------------------
    d/dt(depot) <- transit(nn, mtt, fdepot) - ka * depot
    d/dt(central) <- ka * depot -
      (cl_renal + cl_form_d3og + cl_nonren) * central / vc -
      q * central / vc + q * peripheral1 / vp
    d/dt(peripheral1) <- q * central / vc - q * peripheral1 / vp
    d/dt(central_d3og) <- cl_form_d3og * central / vc -
      cl_d3og * central_d3og / vc_d3og

    # Bioavailability is delivered via transit()'s bio argument
    # above; the depot bolus is suppressed because the entire dose
    # enters the absorption process through the analytical
    # transit-chain input rate (Wilkins 2008 / Tikiso 2021 pattern).
    f(depot) <- 0

    # ------------------------------------------------------------
    # Plasma concentrations (ng/mL with dose in mg and volumes in L
    # the natural rxode2 units yield mg/L which equals ug/mL; the
    # paper-reported scale of ng/mL is obtained by multiplying by
    # 1000). The D3OG state carries dapagliflozin mass-equivalents
    # (1:1 molar from CLP_M15 flux) and V2M / CLM absorb the implicit
    # MW conversion factor per the paper's NONMEM convention -- so
    # the same 1000-x multiplier yields D3OG-mass concentrations in
    # ng/mL that match the paper's plasma D3OG observations.
    # ------------------------------------------------------------
    Cc <- 1000 * central / vc
    Cc_d3og <- 1000 * central_d3og / vc_d3og

    Cc ~ add(addSd) + prop(propSd)
    Cc_d3og ~ add(addSd_d3og) + prop(propSd_d3og)
  })
}
