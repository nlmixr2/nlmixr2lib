Schmid_2017_nintedanib <- function() {
  description <- "Population pharmacokinetic model of nintedanib and its main hydrolytic metabolite BIBF 1202 (Schmid 2017): a one-compartment first-order absorption + lag parent (nintedanib) jointly fit with a one-compartment first-order absorption + lag metabolite (BIBF 1202) coupled to the parent via a fixed fractional formation-during-elimination term (kmet = CL/V2 * ffM) and a fixed V3/V2 volume ratio inherited from rat IV data. The 1191-patient pooled data set spans four trials in NSCLC (Reck 2011 NSCLC phase II, LUME-Lung 1, LUME-Lung 2) and IPF (TOMORROW). Covariates include allometric body weight on CL, linear age on F1, smoking-status on F1, ethnic-origin composite (Indian/Chinese/Taiwanese vs Korean vs reference) on F1, study-group effects on F1 and ka, and on the metabolite side body weight on F2 with ethnic-origin (Indian alone, non-Indian Asian) on F2, ECOG status, LDH (hockey-stick), study-group effect on ka2, and NSCLC histology on ka2."
  reference <- paste(
    "Schmid U, Liesenfeld KH, Fleury A, Dallinger C, Freiwald M.",
    "Population pharmacokinetics of nintedanib, an inhibitor of tyrosine",
    "kinases, in patients with non-small cell lung cancer or idiopathic",
    "pulmonary fibrosis. Cancer Chemotherapy and Pharmacology.",
    "2018 Jan;81(1):89-101. doi:10.1007/s00280-017-3452-0. PMID 29127500.",
    sep = " "
  )
  vignette <- "Schmid_2017_nintedanib"
  units <- list(time = "h", dosing = "mg", concentration = "nmol/L")

  covariateData <- list(
    AGE = list(
      description        = "Age. Linear effect on nintedanib F1 (centered at 62 years).",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline. Reference 62 years (cohort median, Schmid 2017 Table 2). Effect: F1 = ... * (1 + e_age_fdepot * (AGE - 62)) per Schmid 2017 Table 3 footer.",
      source_name        = "AGE"
    ),
    WT = list(
      description        = "Body weight. Power allometric on nintedanib CL and power on BIBF 1202 F2.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline. Reference 71.5 kg (cohort median, Schmid 2017 Table 2). Parent CL: CL = theta_CL * (WT/71.5)^e_wt_cl with e_wt_cl = 0.619 estimated (Schmid 2017 Table 3). BIBF 1202 F2 (relative to F1): F2 = F1 * theta_F2F1 * (WT/71.5)^e_wt_fdepot_bibf with e_wt_fdepot_bibf = -0.848 estimated (Schmid 2017 Online Resource Table S5).",
      source_name        = "WT"
    ),
    SMOKE = list(
      description        = "Current-smoker binary indicator. Multiplicative effect on nintedanib F1 only (no separate effect on BIBF 1202).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ex- or never-smoker; Schmid 2017 cohort: 86.5% of 1191 patients).",
      notes              = "Time-fixed baseline. Effect: F1 = ... * (1 - 0.206)^SMOKE = ... * 0.794^SMOKE, so current smokers have 20.6% lower F1 than ex- or never-smokers (Schmid 2017 Table 3 row 'Current smoker' theta_Smok = 0.794, RSE 4.46%). Schmid 2017 Table 2 reports 'Ex-smoker' (688, 57.8%) and 'Non-smoker' (327, 27.5%) pooled into the SMOKE=0 reference; 'Current smoker' (176, 14.8%) is SMOKE=1.",
      source_name        = "SMOK"
    ),
    RACE_IND_CHI_TWN = list(
      description        = "Indian/Chinese/Taiwanese composite race indicator (1 = Indian, Chinese, or Taiwanese; 0 = Caucasian, Black, Korean, or other Asian). Drives the Indian-Chinese-Taiwanese ethnic-origin effect on nintedanib F1.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Caucasian, Black, Korean, or other Asian; the Schmid 2017 cohort assigns 'other Asian' (subjects of Asian heritage outside China / India / Korea / Taiwan, 3.9% of population) to the reference category alongside Caucasian (75.5%) and Black (0.8%)).",
      notes              = "Time-fixed baseline. Effect on nintedanib F1: F1 = ... * theta_Ethnicity, where theta_Ethnicity = 1.33 if RACE_IND_CHI_TWN = 1, paired exclusively with the RACE_KOREAN = 1 alternative (0.781). Reference 1.00 (Schmid 2017 Table 3). Per the canonical-register entry, this composite is paper-specific to Schmid 2017 and pairs only with the Korean alternative; subjects who are both Indian and (per the operator-distinct RACE_INDIAN canonical) South Asian Indian have both RACE_IND_CHI_TWN = 1 and RACE_INDIAN = 1 -- there is no conflict because the two canonicals address different model-equation slots (RACE_IND_CHI_TWN drives nintedanib F1, RACE_INDIAN drives BIBF 1202 F2).",
      source_name        = "ETHNIC (paper-specific composite of paper-categorical levels Indian / Chinese / Taiwanese)"
    ),
    RACE_KOREAN = list(
      description        = "Korean-heritage race indicator (1 = Korean, 0 = otherwise).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Caucasian, Black, Chinese, Taiwanese, Indian, or other Asian; the Schmid 2017 cohort reference group for the Korean F1 alternative).",
      notes              = "Time-fixed baseline. Effect on nintedanib F1: F1 = ... * 0.781^RACE_KOREAN. Schmid 2017 Table 3 ethnic-origin row Korean theta_Ethnicity = 0.781 (RSE 6.53%). The 5.8% of cohort identified as Korean had ~22% lower F1 than the Caucasian / Black / other-Asian reference. Korean ethnicity is NOT included in the RACE_IND_CHI_TWN composite -- the two canonicals are mutually exclusive: a Korean subject has RACE_IND_CHI_TWN = 0 and RACE_KOREAN = 1. Korean ethnicity ALSO contributes to the RACE_ASIAN broader Asian indicator (= 1 for Korean subjects), which in turn drives the non-Indian-Asian F2 effect through the derived race-asian-non-indian flag inside model().",
      source_name        = "ETHNIC == 'Korean'"
    ),
    RACE_INDIAN = list(
      description        = "South Asian Indian race indicator (1 = Indian, 0 = otherwise). NOT 'American Indian' / Native American.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Caucasian, Black, Chinese, Taiwanese, Korean, or other Asian; the Schmid 2017 cohort reference group for the BIBF 1202 F2 Indian-alone effect).",
      notes              = "Time-fixed baseline. Effect on BIBF 1202 F2: F2 = ... * 1.90^RACE_INDIAN with the non-Indian-Asian alternative 1.20^(RACE_ASIAN AND NOT RACE_INDIAN). Schmid 2017 Online Resource Table S5 theta_Ethnicity for BIBF 1202: 1.90 (RSE 13.5%, Indian alone), 1.20 (RSE 4.83%, non-Indian Asian). The same Indian subjects also have RACE_IND_CHI_TWN = 1 (because the F1 effect groups Indian with Chinese and Taiwanese, a different composite). The two race canonicals do not conflict because they drive distinct model-equation slots (F1 vs F2).",
      source_name        = "ETHNIC == 'Indian'"
    ),
    RACE_ASIAN = list(
      description        = "Broad Asian race indicator (1 = Asian heritage of any subgroup -- Chinese, Korean, Taiwanese, Indian, or other Asian; 0 = non-Asian -- Caucasian or Black). Used in combination with RACE_INDIAN to derive the non-Indian-Asian flag for the BIBF 1202 F2 effect.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Caucasian or Black; the Schmid 2017 cohort reference for BIBF 1202 F2).",
      notes              = "Time-fixed baseline. RACE_ASIAN = 1 for any subject whose RACE_IND_CHI_TWN, RACE_KOREAN, or 'other Asian' bin = 1. Indian subjects also have RACE_ASIAN = 1. The derived flag race_asian_nonind = RACE_ASIAN * (1 - RACE_INDIAN) is computed inside model() and drives the BIBF 1202 F2 effect of 1.20 (RSE 4.83%) per Schmid 2017 Online Resource Table S5.",
      source_name        = "ETHNIC %in% c('Chinese', 'Korean', 'Taiwanese', 'Indian', 'other Asian')"
    ),
    ECOG_GE1 = list(
      description        = "ECOG performance-status indicator (1 = ECOG >= 1, 0 = ECOG = 0). NSCLC patients only; IPF patients have no ECOG and are assigned ECOG_GE1 = 0 by convention (Schmid 2017 model treats IPF subjects as having the reference ECOG status). Effect on BIBF 1202 F2 only (no effect on nintedanib F1).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ECOG = 0; ~22.6% of Schmid 2017 cohort).",
      notes              = "Time-fixed baseline. Effect on BIBF 1202 F2: F2 = ... * 1.16^ECOG_GE1 per Schmid 2017 Online Resource Table S5 (theta_ECOG = 1.16, RSE 4.16%). IPF patients have ECOG missing (Schmid 2017 Table 2 'Missing (due to IPF indication)' 342, 28.7%); in the source model these are pooled with ECOG = 0 reference.",
      source_name        = "ECOG (ECOG_GE1 = as.integer(ECOG >= 1); for IPF cohort with ECOG missing, ECOG_GE1 = 0).",
      reference          = 0
    ),
    LDH = list(
      description        = "Serum lactate dehydrogenase. Hockey-stick effect on BIBF 1202 F2 with estimated breakpoint at 688 U/L; flat above the breakpoint.",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline (Schmid 2017 used LDH at start of treatment). Effect on BIBF 1202 F2: F2 = ... * (1 - e_ldh_coef_fdepot_bibf * (e_ldh_bp_fdepot_bibf - LDH)) when LDH < e_ldh_bp_fdepot_bibf, else F2 = ... * 1. Schmid 2017 Online Resource Table S5 estimates: e_ldh_coef_fdepot_bibf = 0.000656 (RSE 22.1%) and e_ldh_bp_fdepot_bibf = 688 U/L (RSE 22.0%). Cohort median 238 U/L, 5th-95th 141-576 U/L (Schmid 2017 Table 2).",
      source_name        = "LDH"
    ),
    TUMTP_NSCLC_NONADENO = list(
      description        = "NSCLC non-adenocarcinoma histology indicator (1 = NSCLC with non-adenocarcinoma histology, 0 = NSCLC adenocarcinoma OR IPF OR NSCLC of unknown histology). Effect on BIBF 1202 ka2 only (no effect on parent nintedanib ka).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (adenocarcinoma NSCLC, unknown-histology NSCLC, or IPF; pooled per Schmid 2017 Online Resource Table S5 reference).",
      notes              = "Time-fixed baseline. Effect on BIBF 1202 ka2: ka2 = ... * 1.36^TUMTP_NSCLC_NONADENO per Schmid 2017 Online Resource Table S5 (theta_NSCLC_histology = 1.36, RSE 8.90%). Schmid 2017 Table 2: NSCLC non-adenocarcinoma 274 (23.0%); NSCLC adenocarcinoma 502 (42.1%); IPF or NSCLC unknown histology 415 (34.8%).",
      source_name        = "NSCLC_histology == 'non-adenocarcinoma'"
    ),
    STUDY_TOMORROW = list(
      description        = "TOMORROW (NCT00514683; BI 1199.30; Richeldi 2011 NEJM) IPF phase II study indicator. Effect on nintedanib ka and on BIBF 1202 ka2 (paired with STUDY_NSCLC_NIN_PH2).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-TOMORROW: NSCLC phase II, LUME-Lung 1, or LUME-Lung 2).",
      notes              = "Time-fixed per subject. Combined with STUDY_NSCLC_NIN_PH2 to drive the ka phase-II trial-group effect (theta_Trial = 2.20 on nintedanib ka; theta_Trial = 0.756 on BIBF 1202 ka2 per Schmid 2017 Online Resource Table S5).",
      source_name        = "STUDY == 'TOMORROW' (or per-subject study identifier resolving to the IPF phase II / TOMORROW trial)"
    ),
    STUDY_NSCLC_NIN_PH2 = list(
      description        = "Nintedanib phase II NSCLC study indicator (Reck 2011 Ann Oncol; BI 1199.4). Distinct effects on nintedanib F1 (with STUDY_LUMELUNG2) and on nintedanib ka / BIBF 1202 ka2 (with STUDY_TOMORROW).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-NSCLC-phase-II: TOMORROW, LUME-Lung 1, or LUME-Lung 2).",
      notes              = "Time-fixed per subject. Combined with STUDY_LUMELUNG2 to drive the F1 trial-group effect (theta_Trial = 1.30 on F1 per Schmid 2017 Table 3). Combined with STUDY_TOMORROW to drive the ka phase-II trial-group effect (theta_Trial = 2.20 on nintedanib ka per Schmid 2017 Table 3; theta_Trial = 0.756 on BIBF 1202 ka2 per Schmid 2017 Online Resource Table S5).",
      source_name        = "STUDY == 'NSCLC_PHASE2_NINTEDANIB' (Reck 2011 BI 1199.4 phase II)"
    ),
    STUDY_LUMELUNG2 = list(
      description        = "LUME-Lung 2 (NCT00806819; BI 1199.14; Hanna 2016 Lung Cancer) phase III NSCLC + pemetrexed study indicator. Combined with STUDY_NSCLC_NIN_PH2 to drive the F1 trial-group effect.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-LUME-Lung-2).",
      notes              = "Time-fixed per subject. Combined with STUDY_NSCLC_NIN_PH2 to drive the F1 trial-group effect (theta_Trial = 1.30 on F1 per Schmid 2017 Table 3). For the ka / ka2 trial-effect this study acts as part of the reference cohort (paired with LUME-Lung 1).",
      source_name        = "STUDY == 'LUME-Lung_2'"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 1191L,
    n_studies      = 4L,
    age_range      = "45-76 years (5th-95th percentile; median 62)",
    age_median     = "62 years",
    weight_range   = "50.0-100.0 kg (5th-95th percentile; median 71.5)",
    weight_median  = "71.5 kg",
    sex_female_pct = 30.8,
    race_ethnicity = c(Caucasian = 75.5, Asian = 23.7, Black = 0.8),
    asian_breakdown = c(Chinese = 8.2, Korean = 5.8, Indian = 4.2, Taiwanese = 1.6, OtherAsian = 3.9),
    disease_state  = "Advanced non-small cell lung cancer (NSCLC, 71.3%) or idiopathic pulmonary fibrosis (IPF, 28.7%). NSCLC histology: adenocarcinoma 42.1%, non-adenocarcinoma 23.0%, unknown 34.8% (pooled with IPF).",
    dose_range     = "Oral nintedanib 50-250 mg, once- or twice-daily, in repeated treatment cycles continued until disease progression or intolerable toxicity.",
    smoking_distribution = c(NonSmoker_pct = 27.5, ExSmoker_pct = 57.8, CurrentSmoker_pct = 14.8),
    ecog_distribution = c(ECOG_0_pct = 22.6, ECOG_1_pct = 47.2, ECOG_2_pct = 1.5, Missing_IPF_pct = 28.7),
    notes          = "Pooled population PK analysis combining 4 trials: TOMORROW IPF phase II (Richeldi 2011, n = 342 IPF), Reck 2011 nintedanib NSCLC phase II (n = 73 NSCLC, monotherapy), LUME-Lung 1 (Reck 2014, n = 652 NSCLC + docetaxel 75 mg/m^2 q3w), LUME-Lung 2 (Hanna 2016, n = 347 NSCLC + pemetrexed 500 mg/m^2 q3w). Pharmacokinetic data: 5611 nintedanib and 5376 BIBF 1202 plasma concentrations (HPLC-MS/MS), reported in nM (1 nM nintedanib = 1.853 ng/mL; 1 nM BIBF 1202 = 1.903 ng/mL)."
  )

  # Implementation notes (see vignette 'Assumptions and deviations' for the
  # full justification of each item):
  # * Two-output joint parent + metabolite popPK. Compartments are
  #   `depot` and `central` for nintedanib (canonical) and `depot_bibf`
  #   and `central_bibf` for BIBF 1202 (BIBF 1202 added to the
  #   metabolite-suffix register alongside this extraction). The dose
  #   enters BOTH `depot` and `depot_bibf` in parallel; the BIBF 1202
  #   depot receives a small fraction of the dose (theta_F2F1 ~ 0.011 of
  #   F1) representing first-pass conversion of nintedanib to BIBF 1202
  #   in the intestine / liver. A small additional flux of BIBF 1202 is
  #   formed systemically from circulating nintedanib at rate
  #   kmet * central = (CL/V2 * ffM) * central; ffM is fixed at the
  #   IV-bioavailability-study value (Schmid 2017 Online Resource Table
  #   S4, derived from Dallinger 2016 absolute-bioavailability study).
  # * Parameter values are point estimates from the SEQUENTIAL fit
  #   (Schmid 2017 Table 3 for nintedanib; Online Resource Table S5 for
  #   BIBF 1202). The simultaneous-fit alternative (Online Resource
  #   Table S7) is reported as a consistency check and is not the
  #   primary final model.
  # * Concentrations in nM. Conversion factors are MW-derived
  #   (nintedanib MW = 539.62 g/mol -> nM = 1.853 * ng/mL = 1853 * mg/L;
  #   BIBF 1202 MW = 525.59 g/mol -> nM = 1.903 * ng/mL = 1903 * mg/L).
  # * IIV in ka: Schmid 2017 estimated TWO separate IIV values keyed by
  #   trial phase (32.4% CV in phase II vs 53.8% CV in phase III). The
  #   authors state that the phase II estimate is the proper
  #   characterization because phase III sampling was too sparse; the
  #   packaged model uses the phase II value (omega^2 = log(1 + 0.324^2)
  #   = 0.0999) as the single etalka. Documented in vignette Errata.
  # * IIV in V2: 119% CV (Schmid 2017 Table 3) is large but is what the
  #   model estimated. The eta-shrinkage on V2 was 60.0%; we propagate
  #   this IIV value as-is.
  # * Trial-effect covariates: Schmid 2017 encodes between-trial
  #   differences in nintedanib F1 (NSCLC phase II + LUME-Lung 2 vs
  #   IPF phase II + LUME-Lung 1), in nintedanib ka (NSCLC phase II +
  #   IPF phase II vs LUME-Lung 1 + LUME-Lung 2), and in BIBF 1202 ka2
  #   (NSCLC phase II + IPF phase II vs LUME-Lung 1 + LUME-Lung 2;
  #   reverse direction). These are encoded as derived flags inside
  #   model() built from the four mutually exclusive STUDY_* indicators.
  # * Ethnic-origin coverage: Schmid 2017 uses a paper-specific
  #   composite (Indian / Chinese / Taiwanese) on F1 paired with Korean
  #   alone; on F2 the composite is different (Indian alone paired with
  #   non-Indian Asian). The RACE_IND_CHI_TWN canonical handles the F1
  #   composite; RACE_INDIAN handles the F2 Indian-only effect; the
  #   non-Indian Asian F2 flag is derived inside model() from
  #   RACE_ASIAN AND NOT RACE_INDIAN. Operator-ratified sidecar
  #   2026-06-21 (request-002 Q1, option B) instructed full-fidelity
  #   extraction with new race canonicals.
  # * BIBF 1202 ka2 = ka * theta_ka2ka * theta_Trial_ka2 *
  #   theta_NSCLC_histology, so the metabolite ka inherits the parent
  #   ka covariate effects (including the phase-II ka boost = 2.20x)
  #   along with its own additional trial multiplier (0.756x) and
  #   non-adenocarcinoma histology multiplier (1.36x). The model code
  #   reflects this nesting by computing ka first (with parent
  #   covariates), then ka_bibf as a multiplicative deviation.
  # * BIBF 1202 V3/V2 fixed at 0.0185 (Schmid 2017 Online Resource Table
  #   S5; from rat data per Schmid 2017 Online Resource Table S4 and
  #   reference [14]). BIBF 1202 ALAG2 fixed at the parent ALAG (Online
  #   Resource Table S5 theta_ALAG2ALAG = 1.00 fixed). ffM fixed at
  #   0.000931 (Online Resource Table S5; derived from healthy-volunteer
  #   IV nintedanib data per Online Resource Table S4).
  # * Residual error: Schmid 2017 used additive error on the
  #   log-transformed concentration (NONMEM "additive on log scale"),
  #   which in nlmixr2 maps to the lognormal lnorm() residual with the
  #   canonical `expSd` parameter for the parent and `expSd_bibf` for
  #   the metabolite.
  ini({
    # ----- Nintedanib structural parameters (Schmid 2017 Table 3) -----
    lka      <- log(0.0376); label("Log nintedanib absorption rate constant ka (1/h)")    # Schmid 2017 Table 3: theta_ka = 0.0376 (RSE 7.77%; 95% CI 0.0323-0.0439); reference is LUME-Lung 1 + LUME-Lung 2 phase III cohorts (theta_Trial = 1.00)
    lcl      <- log(897);    label("Log nintedanib apparent clearance CL/F at 71.5 kg (L/h)")  # Schmid 2017 Table 3: theta_CL = 897 (RSE 2.42%; 95% CI 855-941)
    lvc      <- log(465);    label("Log nintedanib apparent central volume V2/F (L)")     # Schmid 2017 Table 3: theta_V2 = 465 (RSE 10.7%; 95% CI 376-569)
    ltlag    <- log(0.417);  label("Log nintedanib absorption lag time ALAG (h)")         # Schmid 2017 Table 3: ALAG = 0.417 (RSE 5.59%; 95% CI 0.351-0.463)
    lfdepot  <- fixed(log(1)); label("Log nintedanib reference bioavailability F1 (unitless; fixed = 1)")  # Schmid 2017 Table 3 footnote b: F1 reference fixed = 1

    # ----- Nintedanib covariate effects (Schmid 2017 Table 3) -----
    e_wt_cl                <- 0.619;     label("Power exponent of WT/71.5 on nintedanib CL (unitless)")            # Schmid 2017 Table 3: theta_WT = 0.619 (RSE 16.5%; 95% CI 0.453-0.789)
    e_age_fdepot           <- 0.00959;   label("Linear coefficient on (AGE - 62) for nintedanib F1 (1/year)")     # Schmid 2017 Table 3: theta_Age = 0.00959 (RSE 16.0%; 95% CI 0.00635-0.0126)
    e_smoke_fdepot         <- 0.794;     label("Multiplicative effect on nintedanib F1 for current smokers")      # Schmid 2017 Table 3: theta_Smok 'Current smoker' = 0.794 (RSE 4.46%; 95% CI 0.725-0.864); reference 'Ex- or non-smoker' = 1.00 (fixed)
    e_indchitwn_fdepot     <- 1.33;      label("Multiplicative effect on nintedanib F1 for Indian/Chinese/Taiwanese")  # Schmid 2017 Table 3: theta_Ethnicity 'Indian/Chinese/Taiwanese' = 1.33 (RSE 5.21%; 95% CI 1.19-1.47); reference 'Caucasian/Black/other Asian' = 1.00 (fixed)
    e_korean_fdepot        <- 0.781;     label("Multiplicative effect on nintedanib F1 for Korean origin")        # Schmid 2017 Table 3: theta_Ethnicity 'Korean' = 0.781 (RSE 6.53%; 95% CI 0.690-0.893)
    e_studyhi_fdepot       <- 1.30;      label("Multiplicative effect on nintedanib F1 for NSCLC phase II / LUME-Lung 2 trials")  # Schmid 2017 Table 3: theta_Trial F1 'NSCLC Phase II [23] and LUME-Lung 2 [6]' = 1.30 (RSE 3.77%; 95% CI 1.21-1.39); reference 'IPF Phase II [11] and LUME-Lung 1 [5]' = 1.00 (fixed)
    e_studyphii_ka         <- 2.20;      label("Multiplicative effect on nintedanib ka for NSCLC phase II / IPF phase II trials")  # Schmid 2017 Table 3: theta_Trial ka 'NSCLC Phase II [23] and IPF Phase II [11]' = 2.20 (RSE 8.00%; 95% CI 1.87-2.58); reference 'LUME-Lung 1 [5] and LUME-Lung 2 [6]' = 1.00 (fixed)

    # ----- BIBF 1202 structural parameters (Schmid 2017 Online Resource Table S5) -----
    # Encoded as absolute metabolite parameters derived from the paper's
    # ratio parameterization (kept in comments for source-trace):
    #   theta_F2F1 = 0.0110 -> F2_reference_subject = F1_ref * 0.0110 = 0.0110
    #   theta_ka2ka = 1.47  -> ka_bibf_reference = ka_ref * 1.47 = 0.0553 (1/h)
    #   theta_V3V2 = 0.0185 -> vc_bibf_reference = vc_ref * 0.0185 = 8.6025 L (FIXED)
    #   theta_ALAG2ALAG = 1.00 (FIXED) -> tlag_bibf = tlag (FIXED ratio)
    #   theta_CL2 = 7.05 L/h -> BIBF 1202 apparent CL2/F (estimated)
    #   theta_ffM = 0.000931 (FIXED from healthy-volunteer IV PK, Online Resource Table S4)
    lcl_bibf  <- log(7.05);     label("Log BIBF 1202 apparent clearance CL2/F (L/h)")            # Schmid 2017 Online Resource Table S5: theta_CL2 = 7.05 (RSE 6.37%; 95% CI 6.17-7.93)
    lvc_bibf  <- fixed(log(8.6025));  label("Log BIBF 1202 apparent central volume V3/F at reference (L; fixed via V3/V2 = 0.0185)")  # Schmid 2017 Online Resource Table S5: theta_V3V2 = 0.0185 (fixed; from rat data ref [14]); V3 = V2 * 0.0185 = 465 * 0.0185 = 8.6025
    lka_bibf  <- log(0.0553);   label("Log BIBF 1202 absorption rate constant ka2 reference (1/h)")  # Schmid 2017 Online Resource Table S5: theta_ka2ka = 1.47 (RSE 8.16%; 95% CI 1.23-1.71); ka2 = ka * 1.47 = 0.0376 * 1.47 = 0.0553 (1/h)
    ltlag_bibf <- fixed(log(0.417));  label("Log BIBF 1202 absorption lag time ALAG2 (h; fixed = parent ALAG)")  # Schmid 2017 Online Resource Table S5: theta_ALAG2ALAG = 1.00 (fixed; ALAG2 = ALAG = 0.417 h)
    lfdepot_bibf <- log(0.0110);  label("Log BIBF 1202 reference bioavailability F2 (unitless; = F1 * 0.0110)")  # Schmid 2017 Online Resource Table S5: theta_F2F1 = 0.0110 (RSE 11.3%; 95% CI 0.00856-0.0134); F2 = F1 * 0.0110 at reference subject
    lffm      <- fixed(log(0.000931));  label("Log fraction of nintedanib systemic CL forming BIBF 1202 (unitless; fixed)")  # Schmid 2017 Online Resource Table S5: theta_ffM = 0.000931 (fixed); derived from healthy-volunteer IV nintedanib PK (Online Resource Table S4)

    # ----- BIBF 1202 covariate effects (Schmid 2017 Online Resource Table S5) -----
    e_wt_fdepot_bibf       <- -0.848;    label("Power exponent of WT/71.5 on BIBF 1202 F2 (unitless)")             # Schmid 2017 Online Resource Table S5: theta_WT (F2) = -0.848 (RSE 12.2%; 95% CI -0.646 to -1.05)
    e_indian_fdepot_bibf   <- 1.90;      label("Multiplicative effect on BIBF 1202 F2 for Indian origin alone")   # Schmid 2017 Online Resource Table S5: theta_Ethnicity 'Indian' = 1.90 (RSE 13.5%; 95% CI 1.40-2.40); reference 'Caucasian/Black' = 1.00 (fixed)
    e_asiannonind_fdepot_bibf <- 1.20;   label("Multiplicative effect on BIBF 1202 F2 for non-Indian Asian origin")  # Schmid 2017 Online Resource Table S5: theta_Ethnicity 'Asian except Indian' = 1.20 (RSE 4.83%; 95% CI 1.09-1.31)
    e_ldh_coef_fdepot_bibf <- 0.000656;  label("Hockey-stick LDH slope on BIBF 1202 F2 (1/(U/L))")                # Schmid 2017 Online Resource Table S5: theta_LDH = 0.000656 (RSE 22.1%; 95% CI 0.000372-0.000940)
    e_ldh_bp_fdepot_bibf   <- 688;       label("Hockey-stick LDH breakpoint on BIBF 1202 F2 (U/L)")               # Schmid 2017 Online Resource Table S5: theta_LDH_breakpoint = 688 (RSE 22.0%; 95% CI 391-985)
    e_ecog_fdepot_bibf     <- 1.16;      label("Multiplicative effect on BIBF 1202 F2 for ECOG >= 1")             # Schmid 2017 Online Resource Table S5: theta_ECOG 'ECOG >= 1' = 1.16 (RSE 4.16%; 95% CI 1.07-1.25); reference 'ECOG = 0' = 1.00 (fixed)
    e_studyphii_ka_bibf    <- 0.756;     label("Multiplicative effect on BIBF 1202 ka2 for NSCLC phase II / IPF phase II trials")  # Schmid 2017 Online Resource Table S5: theta_Trial (ka2) = 0.756 (RSE 7.72%; 95% CI 0.642-0.870); reference 'LUME-Lung 1 + LUME-Lung 2' = 1.00 (fixed)
    e_nonaden_ka_bibf      <- 1.36;      label("Multiplicative effect on BIBF 1202 ka2 for non-adenocarcinoma NSCLC")  # Schmid 2017 Online Resource Table S5: theta_NSCLC_histology 'non-adenocarcinoma' = 1.36 (RSE 8.90%; 95% CI 1.12-1.60); reference 'adenocarcinoma' = 1.00 (fixed)

    # ----- Inter-individual variability (Schmid 2017 Table 3 and Online Resource Table S5) -----
    # Convert reported CV% to omega^2 via the exact log-normal identity
    # omega^2 = log(1 + CV^2). Schmid 2017 reports CV in percent on the
    # original scale; the conversion preserves the geometric CV when
    # the model uses log-normal parameterisation.
    etalfdepot       ~ 0.21500   # Schmid 2017 Table 3: IIV in F1 = 49.1% CV; log(1 + 0.491^2) = 0.21500
    etalka           ~ 0.09989   # Schmid 2017 Table 3: IIV in ka (phase II) = 32.4% CV; log(1 + 0.324^2) = 0.09989. Phase III value (53.8% CV) is NOT used; the authors state that the phase II estimate is the proper characterization (Schmid 2017 Table 3 footnote d).
    etalvc           ~ 0.86340   # Schmid 2017 Table 3: IIV in V2/F = 119% CV; log(1 + 1.19^2) = 0.86340
    etalfdepot_bibf  ~ 0.28045   # Schmid 2017 Online Resource Table S5: IIV in F2 = 56.7% CV; log(1 + 0.567^2) = 0.28045

    # ----- Residual unexplained variability (Schmid 2017 Table 3 and Online Resource Table S5) -----
    # Schmid 2017 reports residual error as additive on the
    # log-transformed concentration scale (NONMEM "additive on log");
    # this maps in nlmixr2 to the lognormal residual via lnorm(expSd)
    # with expSd being the SD on the natural-log scale.
    expSd        <- 0.526; label("Lognormal residual SD on nintedanib plasma concentration (log[nM])")  # Schmid 2017 Table 3: additive SD = 0.526 [nM; log scale] (RSE 4.58%; 95% CI 0.504-0.553)
    expSd_bibf   <- 0.546; label("Lognormal residual SD on BIBF 1202 plasma concentration (log[nM])")  # Schmid 2017 Online Resource Table S5: additive SD = 0.546 [nM; log scale] (RSE 4.13%; 95% CI 0.523-0.568)
  })
  model({
    # ----- Reference covariate values (Schmid 2017 Table 2 medians) -----
    ref_wt  <- 71.5
    ref_age <- 62

    # ----- Concentration unit conversion factors -----
    # central in mg, vc in L -> mg/L = 1000 ng/mL. Schmid 2017 reports
    # plasma concentrations in nM with conversion factors derived from
    # the molecular weights of nintedanib (MW = 539.62) and BIBF 1202
    # (MW = 525.59): 1 nM nintedanib = 1.853 ng/mL; 1 nM BIBF 1202 =
    # 1.903 ng/mL (Schmid 2017 Table 3 footnote and Online Resource
    # Table S7 cont. footnote).
    cf_nint <- 1853       # nM per mg/L for nintedanib
    cf_bibf <- 1903       # nM per mg/L for BIBF 1202

    # ----- Derived trial-group flags (mutually exclusive STUDY_* indicators) -----
    # A subject is in exactly one of the four studies; the OR-sums are
    # therefore 0 or 1 (never 2). The reference categories below
    # implicitly assign 'all four indicators = 0' to the LUME-Lung 1 +
    # something else; the model uses the explicit OR-sums to be robust
    # to coding accidents.
    trial_F1_hi  <- STUDY_NSCLC_NIN_PH2 + STUDY_LUMELUNG2   # F1 trial-group high (NSCLC phase II OR LUME-Lung 2)
    trial_KA_ph2 <- STUDY_NSCLC_NIN_PH2 + STUDY_TOMORROW    # ka phase-II trial-group (NSCLC phase II OR IPF phase II)

    # ----- Derived race flag for non-Indian Asian (BIBF 1202 F2) -----
    race_asian_nonind <- RACE_ASIAN * (1 - RACE_INDIAN)

    # ----- Individual nintedanib parameters -----
    # Allometric WT on CL (estimated exponent 0.619).
    cl <- exp(lcl) * (WT / ref_wt)^e_wt_cl

    # V2 with IIV (no covariates other than the log-normal random effect).
    vc <- exp(lvc + etalvc)

    # Absorption rate ka with IIV and the phase-II trial multiplier.
    ka <- exp(lka + etalka) * e_studyphii_ka^trial_KA_ph2

    # Absorption lag time (no IIV, no covariates).
    tlag <- exp(ltlag)

    # Nintedanib relative bioavailability F1 with covariates:
    #   - Ethnicity (Indian/Chinese/Taiwanese, Korean)
    #   - Age (linear in AGE - 62)
    #   - Smoking (current smoker)
    #   - Trial group (NSCLC phase II + LUME-Lung 2 vs reference)
    #   - IIV on F1
    fdepot <- exp(lfdepot + etalfdepot) *
              e_indchitwn_fdepot^RACE_IND_CHI_TWN *
              e_korean_fdepot^RACE_KOREAN *
              (1 + e_age_fdepot * (AGE - ref_age)) *
              e_smoke_fdepot^SMOKE *
              e_studyhi_fdepot^trial_F1_hi

    # ----- Individual BIBF 1202 parameters -----
    # CL2/F estimated, no IIV reported.
    cl_bibf <- exp(lcl_bibf)

    # V3/F = V2/F * theta_V3V2 (theta_V3V2 = 0.0185 fixed). The model
    # encodes lvc_bibf = log(V3 at reference) = log(8.6025) FIXED, so
    # V3 inherits IIV from V2 via the same etalvc (the ratio is fixed,
    # so individual deviations in V2 propagate proportionally to V3).
    vc_bibf <- exp(lvc_bibf + etalvc)

    # ka2 = ka * theta_ka2ka * theta_Trial_ka2 * theta_NSCLC_histology
    # The parent ka already carries etalka and the phase-II trial
    # multiplier (2.20); ka2 adds its own multipliers on top:
    #   - theta_ka2ka = 1.47 (encoded as exp(lka_bibf - lka))
    #   - theta_Trial_ka2 = 0.756 for NSCLC phase II + IPF phase II
    #   - theta_NSCLC_histology = 1.36 for non-adenocarcinoma NSCLC
    # The model body computes ka2 by multiplying the parent ka by the
    # encoded ratio and applying the metabolite-specific covariates.
    ka_bibf <- ka * exp(lka_bibf - lka) *
               e_studyphii_ka_bibf^trial_KA_ph2 *
               e_nonaden_ka_bibf^TUMTP_NSCLC_NONADENO

    # ALAG2 = ALAG * theta_ALAG2ALAG (theta_ALAG2ALAG = 1.00 fixed); no
    # additional covariates. The encoded ltlag_bibf is fixed at the
    # parent ltlag value, so tlag_bibf is identically tlag.
    tlag_bibf <- exp(ltlag_bibf)

    # BIBF 1202 hockey-stick LDH factor: (1 - coef * (bp - LDH)) when
    # LDH < bp, else 1. The (LDH < e_ldh_bp_fdepot_bibf) term returns
    # 1 if true and 0 if false in rxode2, so multiplying by that flag
    # zeros out the correction when LDH is at or above the breakpoint.
    ldh_below_bp <- (LDH < e_ldh_bp_fdepot_bibf)
    ldh_drop     <- (e_ldh_bp_fdepot_bibf - LDH) * ldh_below_bp
    ldh_factor   <- 1 - e_ldh_coef_fdepot_bibf * ldh_drop

    # BIBF 1202 F2 = F1 * theta_F2F1 * (WT/71.5)^theta_WT *
    #   theta_Ethnicity_BIBF * theta_ECOG * ldh_factor * exp(eta_F2)
    # The model uses the encoded lfdepot_bibf = log(theta_F2F1)
    # together with the etalfdepot_bibf IIV on top of the F1
    # (parent F1 already includes etalfdepot via the fdepot factor).
    fdepot_bibf <- fdepot *
                   exp(lfdepot_bibf + etalfdepot_bibf) *
                   (WT / ref_wt)^e_wt_fdepot_bibf *
                   e_indian_fdepot_bibf^RACE_INDIAN *
                   e_asiannonind_fdepot_bibf^race_asian_nonind *
                   e_ecog_fdepot_bibf^ECOG_GE1 *
                   ldh_factor

    # ----- Metabolite formation rate constant -----
    # kmet = (CL/V2) * ffM. Schmid 2017: ffM is the small fraction
    # (0.000931, fixed) of nintedanib elimination that produces
    # systemic BIBF 1202. The bulk of BIBF 1202 enters via the
    # dose-routed depot_bibf compartment (first-pass formation,
    # captured by F2 ~ 0.011 of F1); kmet adds a slow systemic
    # formation flux from the central nintedanib compartment.
    ffm  <- exp(lffm)
    kmet <- (cl / vc) * ffm

    # ----- Bioavailability assignments -----
    f(depot)      <- fdepot
    f(depot_bibf) <- fdepot_bibf

    # ----- Absorption lag-time assignments -----
    alag(depot)      <- tlag
    alag(depot_bibf) <- tlag_bibf

    # ----- ODE system -----
    # Parent nintedanib: first-order absorption from depot into central,
    # linear elimination from central.
    d/dt(depot)        <- -ka * depot
    d/dt(central)      <-  ka * depot - (cl / vc) * central

    # BIBF 1202 metabolite: first-order absorption from depot_bibf into
    # central_bibf (capturing first-pass formation); systemic formation
    # at rate kmet * central added to the inflow; linear elimination
    # at rate cl_bibf / vc_bibf.
    d/dt(depot_bibf)   <- -ka_bibf * depot_bibf
    d/dt(central_bibf) <-  ka_bibf * depot_bibf +
                           kmet * central -
                           (cl_bibf / vc_bibf) * central_bibf

    # ----- Observations in nM -----
    Cc      <- (central      / vc)      * cf_nint   # nintedanib (nM)
    Cc_bibf <- (central_bibf / vc_bibf) * cf_bibf   # BIBF 1202 (nM)

    # ----- Residual error: additive on log scale (= lognormal) -----
    Cc      ~ lnorm(expSd)
    Cc_bibf ~ lnorm(expSd_bibf)
  })
}
