Lin_2024_TB_multistate <- function() {
  description <- paste(
    "Pharmacometric five-state multistate model for long-term treatment outcomes in adults with drug-resistant",
    "pulmonary tuberculosis treated with bedaquiline on top of a multidrug background regimen (pooled TMC207-C208",
    "and TMC207-C209 phase IIb trials). Five ODE states carry the marginal state-occupancy probability of active TB",
    "(s_activeTb, the initial state), sputum-culture conversion (s_converted), recurrent TB (s_recurrentTb),",
    "study dropout (s_dropout, absorbing) and death (s_death, absorbing). Patients move from active TB to converted",
    "and on to recurrent TB, and may return from recurrent TB to converted; dropout and death compete from every",
    "transient state. The conversion hazard lambda12 follows a symmetric surge function of time (amplitude SA12,",
    "peak time PT12, surge width SW12) that peaks at 11.4 weeks; recurrence (lambda23), re-conversion (lambda32),",
    "dropout (lambda14, lambda24 = lambda34) and death-from-converted (lambda25) are constant, while death from",
    "active or recurrent TB (lambda15 = lambda35) follows a Weibull hazard with increasing risk over time",
    "(shape 1.96). Faster model-derived bacterial clearance through week 2 raises the surge amplitude and moves the",
    "peak earlier; a longer baseline MGIT time-to-positivity (lower bacterial burden) moves the peak earlier and",
    "widens the surge; XDR-TB lowers the surge amplitude by 46%. Recurrence is higher in men and rises with the",
    "model-derived mycobacterial load at the end of the 24-week treatment period (applied only after week 26).",
    "Dropout is higher in the C208 study and in younger patients; death risk falls with higher baseline body",
    "weight. There are no drug-dosing events and bedaquiline exposure does not enter this model directly -- it",
    "reaches the model through the two model-derived covariates MBL_HL_WK2 and MBL_END, which are secondary metrics",
    "of the upstream mycobacterial-load model shipped as modellib('Svensson_2017_bedaquiline'). The source analysis",
    "was fitted in NONMEM with the exact likelihood on the observed categorical state and has no between-subject",
    "variability (OMEGA fixed to 0) and no residual-error model."
  )
  reference <- paste(
    "Lin YJ, Zou Y, Karlsson MO, Svensson EM.",
    "A pharmacometric multistate model for predicting long-term treatment outcomes of patients with pulmonary TB.",
    "J Antimicrob Chemother. 2024;79(10):2561-2569.",
    "doi:10.1093/jac/dkae256. PMCID: PMC11441995.",
    "Open Access under CC BY.",
    "Structural equations from Supplementary Data 'Model equations' (equations 1-9) and the verbatim NONMEM",
    "control stream in Supplementary Data 'NONMEM code' (pages 12-16); parameter values from Supplementary",
    "Table S2 cross-checked against the control stream $THETA block.",
    "The two model-derived covariates are secondary metrics of Svensson EM, Karlsson MO.",
    "J Antimicrob Chemother. 2017;72(12):3398-3405; see modellib('Svensson_2017_bedaquiline').",
    sep = " "
  )
  vignette <- "Lin_2024_TB_multistate"

  # Multistate state-occupancy probabilities. Same family as the
  # `s_alive` / `s_dropout` / `s_death` states of
  # Ibrahim_2023_ibrutinib_competing_risk.R.
  paper_specific_compartments <- c(
    "s_activeTb", "s_converted", "s_recurrentTb", "s_dropout", "s_death"
  )

  units <- list(
    time          = "week",
    dosing        = "n/a (no drug-dosing events; the model propagates state-occupancy probabilities)",
    concentration = "probability (all five states are occupancy probabilities, not drug concentrations)"
  )

  compartmentData <- list(
    s_activeTb    = list(analyte = "probability of being in the active TB state", units = NA_character_, specimen = "not applicable", verified = TRUE),
    s_converted   = list(analyte = "probability of being in the sputum-culture-converted state", units = NA_character_, specimen = "not applicable", verified = TRUE),
    s_recurrentTb = list(analyte = "probability of being in the recurrent TB state", units = NA_character_, specimen = "not applicable", verified = TRUE),
    s_dropout     = list(analyte = "probability of having dropped out of the study", units = NA_character_, specimen = "not applicable", verified = TRUE),
    s_death       = list(analyte = "probability of having died", units = NA_character_, specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    MBL_HL_WK2 = list(
      description        = "Model-derived half-life of mycobacterial-load decline evaluated over the first 2 weeks of treatment (the source's HL2).",
      units              = "week",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject in the final model. Lin 2024 Supplementary 'Predictors investigation' derives the",
        "metric prospectively from data up to week 2 only ('only data until time t were used to compute",
        "model-derived predictors at time t'), and carries the week-2 value forward for the rest of the study, so",
        "a single scalar per subject reproduces the final model exactly. Enters BOTH surge parameters of the",
        "conversion hazard with a single shared coefficient of opposite sign:",
        "SA12 = exp(lsa12 + e_hl2_sa12pt12 * (MBL_HL_WK2 - 0.69443)) and",
        "PT12 = exp(lpt12 - e_hl2_sa12pt12 * (MBL_HL_WK2 - 0.69443)), with e_hl2_sa12pt12 = -0.686145",
        "(NONMEM code $PK, THETA(12); Table S2 'betaHL2 on SA12, PT12' = -0.686). Because the coefficient is",
        "negative, a shorter half-life (faster bacterial clearance) raises the surge amplitude AND moves the peak",
        "earlier, matching Lin 2024 Discussion paragraph 1. Centering value 0.69443 weeks is the cohort median,",
        "hardcoded in the control stream and confirmed by Lin 2024 Figure 3, whose reference individual has",
        "'0.69 weeks of half-life of bacterial clearance at Week 2'. Lin 2024 Figure 6 gives the 5th and 95th",
        "percentiles as 0.33 weeks (high clearance) and 1.1 weeks (low clearance). The intended source is the",
        "half-life of the mycobacterial-load state of the upstream model",
        "modellib('Svensson_2017_bedaquiline'), whose hl_i is already in weeks -- no conversion needed. Must be",
        "strictly positive."
      ),
      source_name        = "HL2"
    ),
    MBL_END = list(
      description        = "Model-derived mycobacterial load at the end of the 24-week bedaquiline treatment period (the source's MMBLend), on the natural scale.",
      units              = "n bacteria per sample inoculum",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject. Enters the recurrence hazard on the NATURAL-LOG scale and only after week 26:",
        "lambda23 = exp(llambda23 + e_sexf_lambda23 * SEXF + e_mblend_lambda23 * (t > 26) *",
        "(log(MBL_END) - log(5.5726e-05))), with e_mblend_lambda23 = 0.0371081 (NONMEM code $DES, THETA(18);",
        "Table S2 'betaMMBLend on lambda23' = 0.0371). The centering value 5.5726e-05 is hardcoded in the control",
        "stream as LOG(0.000055726) and is the cohort median; log10(5.5726e-05) = -4.25, matching the",
        "'-4.3 log10(MMBLend)' reference individual of Lin 2024 Figure 3. NOTE the coefficient is per unit of",
        "NATURAL log, so the per-log10-unit hazard ratio is exp(0.0371081 * ln(10)) = 1.089. Units match the",
        "mycobacterial-load state of the upstream model modellib('Svensson_2017_bedaquiline') (n bacteria per",
        "sample inoculum, typical starting MBL_0 = 2.14e3), so its mbl state at week 24 can be supplied directly.",
        "Must be strictly positive -- log(0) is undefined and would propagate NaN through the whole solve even",
        "before week 26.",
        "TIME-GATE DISCREPANCY: the Lin 2024 main text says this predictor acts 'after 24 weeks' / 'after the",
        "completion of 24 week treatment', but the control stream gates it at IF(WEEK.GT.26). The control stream",
        "is used here (see the vignette Assumptions and deviations section)."
      ),
      source_name        = "MMBLend"
    ),
    TTP_MGIT_BASE = list(
      description        = "Baseline mean time-to-positivity in the mycobacterial growth indicator tube (MGIT) liquid-culture system, as the mean of triplicate pre-treatment sputum samples.",
      units              = "days",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject. Enters BOTH the peak time and the surge width of the conversion hazard with a",
        "single shared coefficient of opposite sign:",
        "PT12 = exp(lpt12 - e_ttp_pt12sw12 * (TTP_MGIT_BASE - 9.06944) / 7) and",
        "SW12 = exp(lsw12 + e_ttp_pt12sw12 * (TTP_MGIT_BASE - 9.06944) / 7), with",
        "e_ttp_pt12sw12 = 0.442668 (NONMEM code $PK, THETA(11); Table S2 'betabasTTP on PT12, SW12' = 0.443).",
        "A longer baseline TTP (lower bacterial burden) therefore moves the conversion peak earlier and widens the",
        "surge, matching Lin 2024 Results 'Predictors of state transitions'.",
        "UNIT AND CENTERING NOTE: the source data column MTTP is in HOURS, centered at 217.6667 h and then divided",
        "by 24 and by 7 to express the deviation in WEEKS. This canonical is in DAYS, so the identical arithmetic",
        "is written here as (TTP_MGIT_BASE - 9.06944) / 7 with 9.06944 days = 217.6667 / 24; that matches the",
        "'9.1 (2.3-42)' days cohort median of Lin 2024 Table 1. The divisor 7 converts the day deviation to weeks,",
        "so the coefficient is per WEEK of baseline TTP.",
        "COHORT-SPECIFIC REFERENCE: the register entry for this canonical records 6.8 days from the Svensson 2017",
        "C208-only cohort (n = 191). Lin 2024 pools C208 + C209 (n = 402) and its median is 9.06944 days; per the",
        "register's own instruction, the cohort-specific reference is documented here rather than overwriting the",
        "canonical."
      ),
      source_name        = "MTTP"
    ),
    DIS_TB_XDR_STRICT = list(
      description        = "Extensively-drug-resistant (XDR) tuberculosis indicator, XDR only: 1 = XDR, 0 = pre-XDR, MDR, drug-susceptible or unclassified.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-XDR: pre-XDR + MDR + drug-susceptible + missing-treated-as-MDR)",
      notes              = paste(
        "Time-fixed per subject. Encoded in the NONMEM control stream as XDR = 0; IF(TBTYPE.EQ.4) XDR = 1, i.e.",
        "the pre-XDR category (95 of 402 patients, 28%) sits in the REFERENCE group -- this is a strictly narrower",
        "dichotomisation than the sibling canonical DIS_TB_XDR, which pools pre-XDR with XDR. Lin 2024 Results",
        "states the contrast explicitly as XDR-TB 'compared with those with non-XDR-TB'.",
        "Enters the conversion surge amplitude only: SA12 = exp(lsa12 + ... + e_xdr_sa12 * DIS_TB_XDR_STRICT) with",
        "e_xdr_sa12 = -0.622792 (NONMEM code $PK, THETA(13); Table S2 'betaXDR-TB on SA12' = -0.623).",
        "Falsification check: exp(-0.622792) = 0.536, i.e. a 46.4% reduction, matching the paper's",
        "'lambda12 decreased by 46% (95% CI, 24%-62%)'.",
        "Per Lin 2024 Discussion, the ~16% of patients with missing drug-resistance profile were assigned to the",
        "MDR (reference) group after testing that they did not differ significantly from MDR."
      ),
      source_name        = "XDR"
    ),
    SEXF = list(
      description        = "Sex: 1 = female, 0 = male.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "Time-fixed per subject. Enters the recurrence hazard as exp(e_sexf_lambda23 * SEXF) with",
        "e_sexf_lambda23 = -0.812999 (NONMEM code $DES, THETA(17); Table S2 'betasex on lambda23' = -0.813), a",
        "hazard ratio of exp(-0.813) = 0.443 for women relative to men -- i.e. men are at higher risk of",
        "recurrence, matching Lin 2024 Results and Discussion ('the risk of recurrence was higher for men').",
        "POLARITY TRACE: the source column is named SEX and missing values are imputed to 0 (IF(SEX2.EQ.-99)",
        "SEX2 = 0). Because Lin 2024 imputes categorical covariates with the MODE and the cohort is 65% male",
        "(Table 1), the imputed value 0 identifies 0 = male and hence 1 = female, so the source column is already",
        "a female indicator and needs no transformation to reach the SEXF canonical. Independently confirmed by",
        "Lin 2024 Figure 3, whose reference individual (all covariates at their typical value, hence SEXF = 0) is",
        "described as 'a 33-year-old male'."
      ),
      source_name        = "SEX"
    ),
    AGE = list(
      description        = "Age at baseline.",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject. Enters all three dropout hazards (lambda14, lambda24 and lambda34, which share a",
        "single covariate exponent) as exp(e_age_lambda1424 * (AGE - 33)) with e_age_lambda1424 = -0.0230092",
        "(NONMEM code $DES, THETA(15); Table S2 'betaage on lambda14/24/34' = -0.0230). Younger patients drop out",
        "more, matching Lin 2024 Results. Centering value 33 years is the cohort median (Lin 2024 Table 1,",
        "'33 (18-68)') and is hardcoded in the control stream both as the median-imputation value for missing age",
        "and as the centering constant."
      ),
      source_name        = "AGE"
    ),
    WT = list(
      description        = "Body weight at baseline.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject: this model uses BASELINE weight only. Lin 2024 screened time-varying body-weight",
        "change as a predictor (Methods, 'Post-baseline time-varying predictors') but did not retain it, so there",
        "is no time-varying weight column to pair with and the plain WT canonical is used rather than WT_BASE",
        "(which is reserved for the Wahlby 2004 baseline/time-varying split).",
        "Enters all three death hazards (lambda15, lambda25 and lambda35, which share a single covariate exponent)",
        "as exp(e_baswt_lambda1525 * (WT - 55)) with e_baswt_lambda1525 = -0.0838138 (NONMEM code $DES,",
        "THETA(16); Table S2 'betabasWT on lambda15/25/35' = -0.0838). Lower baseline weight carries a higher risk",
        "of death, matching Lin 2024 Results. Centering value 55 kg is the cohort median (Lin 2024 Table 1,",
        "'55 (30-113)'). Lin 2024 Discussion cautions that this finding 'should be interpreted with caution",
        "outside the investigated weight range' of 30-113 kg."
      ),
      source_name        = "WT"
    ),
    STUDY_C208 = list(
      description        = "Source-study indicator for the pooled analysis: 1 = TMC207-C208 (randomized, double-blind, placebo-controlled phase IIb), 0 = TMC207-C209 (open-label, single-arm phase IIb).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (TMC207-C209)",
      notes              = paste(
        "Time-fixed per subject. Enters all three dropout hazards (lambda14, lambda24 and lambda34, which share a",
        "single covariate exponent) as exp(e_study_lambda1424 * STUDY_C208) with e_study_lambda1424 = 0.909188",
        "(NONMEM code $DES, THETA(14); Table S2 'betastudy on lambda14/24/34' = 0.910), a hazard ratio of",
        "exp(0.909188) = 2.48 for C208 relative to C209. Matches Lin 2024 Discussion: 'The dropout rate in patients",
        "from the C208 study was higher than in C209'. The reference is C209, confirmed by Lin 2024 Figure 3, whose",
        "reference individual is 'enrolled in the C209 study'. Cohort split: C208 n = 195 (49%), C209 n = 207 (51%)",
        "(Lin 2024 Table 1)."
      ),
      source_name        = "C208"
    )
  )

  # Covariates that Lin 2024 explicitly screened but did NOT retain in the
  # final model. Documented here for provenance; they are deliberately absent
  # from model() and from covariateData.
  covariatesDataExcluded <- list(
    HIV_POS = list(
      description = "Concomitant HIV infection (Y/N).",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as a baseline covariate on every transition; not retained. Lin 2024 Results: 'The presence of cavitation and HIV status were not identified as predictors of any transition hazard between each state.' The Discussion attributes this to 91% of the cohort not living with HIV (36 of 402, 9%, were HIV-positive; Table 1)."
    ),
    DIS_TB_CAVITATION = list(
      description = "Presence of lung cavitation on baseline chest radiograph (Y/N).",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as a baseline covariate on every transition; not retained (Lin 2024 Results). The Discussion attributes this to 96% of the cohort (385 of 402, Table 1) having cavitary disease, leaving little contrast."
    ),
    DIS_TB_PRIORTX = list(
      description = "Previous anti-TB treatment before enrollment in the study (Y/N).",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as a baseline covariate; not retained. Lin 2024 Discussion notes the information 'was indirectly included through the bacterial clearance and MMBL metrics', i.e. it already acts through MBL_HL_WK2 and MBL_END. 186 of 402 patients (46%) had prior treatment (Table 1)."
    ),
    TRT_BDQ = list(
      description = "Received bedaquiline on top of the background regimen (Y/N); 303 of 402 patients (75%).",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as a baseline covariate; not retained as a direct predictor of any transition. Bedaquiline exposure reaches the model indirectly through MBL_HL_WK2 and MBL_END, which are derived from the exposure-driven upstream mycobacterial-load model modellib('Svensson_2017_bedaquiline')."
    ),
    ALB = list(
      description = "Serum albumin over time, derived from the upstream bedaquiline population PK model.",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened as a post-baseline time-varying predictor (Lin 2024 Methods); not retained in the final model. Computed in the source from a self-limiting logistic model in Svensson 2016 (CPT PSP)."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 402L,
    n_studies      = 2L,
    age_range      = "18-68 years (median 33)",
    age_median     = "33 years",
    weight_range   = "30-113 kg (median 55)",
    weight_median  = "55 kg",
    sex_female_pct = 35,
    race_ethnicity = NULL,
    disease_state  = "Drug-resistant pulmonary tuberculosis: 3% drug-sensitive, 54% MDR, 28% pre-XDR, 14% XDR (pre-2021 WHO definitions); 96% with lung cavitation; 9% living with HIV; 46% previously treated for TB",
    dose_range     = "Bedaquiline 400 mg once daily for 2 weeks, then 200 mg three times weekly for a further 22 weeks (6 weeks in C208 stage 1), on top of a multidrug background regimen; 75% of patients received bedaquiline and 25% received placebo (C208 stage 1/2 placebo arm)",
    regions        = "Brazil, China, Estonia, India, Latvia, Peru, Philippines, Russia, South Africa, South Korea, Thailand, Turkey, Ukraine; 75% of patients from a high-TB-burden country",
    notes          = paste(
      "Baseline characteristics from Lin 2024 Table 1. Of 439 patients enrolled in TMC207-C208 (NCT00449644) and",
      "TMC207-C209 (NCT00910871), 35 with negative cultures at both screening and baseline and 2 with only",
      "pre-treatment observations were excluded, leaving 402. Nineteen patients had missing baseline sputum",
      "culture information and were assumed to start in the active TB state. Follow-up was 120 weeks (104 weeks",
      "for C208 stage 1): a 24-week investigational treatment period plus a 96-week follow-up period.",
      "The dataset holds 6984 state observations with 364 conversion events (25 of them a second conversion after",
      "recurrence), 63 recurrence events, 114 dropouts and 30 deaths (Lin 2024 Figure 1 and Results).",
      "The multistate model was fitted in NONMEM 7.5.1 with the exact likelihood method and $OMEGA 0 FIX, so it",
      "carries no between-subject variability and no residual-error model."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # All hazards are on a WEEK time axis (units$time above).
    #
    # SOURCE-TRACE NOTE ON UNITS. The Lin 2024 NONMEM control stream runs on
    # an HOUR time axis and parameterises the constant hazards as mean
    # transit times in weeks that are converted to hours (THETA(n)*24*7).
    # Supplementary Table S2 reports the same quantities as RATES in
    # week^-1. This file uses Table S2's week^-1 rate parameterisation, and
    # every comment below gives BOTH the Table S2 rounded rate and the
    # unrounded control-stream $THETA it is derived from, so the arithmetic
    # is auditable. Worked example: THETA(1) = 295.375 weeks, so
    # lambda23 = 1/295.375 = 0.00338554 week^-1, which rounds to the
    # 0.00339 printed in Table S2.
    # ------------------------------------------------------------------

    # --- Conversion hazard lambda12: symmetric surge function ----------
    # Supplement equation 8: lambda12(t) = SA / (((t - PT)/SW)^2 + 1)
    lsa12 <- log(0.18967032)
    label("Log surge amplitude SA12 of the conversion hazard (1/week)")  # Table S2 SA12 = 0.190 /week (RSE 9.5%); = THETA(7)/10000*24*7 = 11.2899/10000 * 168 = 0.18967032 /week (control stream states SA in hours^-1)
    lpt12 <- log(11.3799)
    label("Log peak time PT12 of the conversion surge (week)")  # Table S2 PT12 = 11.4 weeks (RSE 6.4%); $THETA 8 = 11.3799. Paper Results: peak at 11.4 weeks (95% CI 10.0-12.8)
    lsw12 <- log(5.58883)
    label("Log surge width SW12 of the conversion surge (week)")  # Table S2 SW12 = 5.59 weeks (RSE 12%); $THETA 9 = 5.58883

    # --- Constant transition hazards -----------------------------------
    llambda23 <- log(0.0033855427)
    label("Log baseline hazard of recurrence, converted -> recurrent TB (1/week)")  # Table S2 lambda23 = 0.00339 /week (RSE 14%); = 1/295.375 ($THETA 1 MTT23) = 0.0033855427
    llambda32 <- log(0.016162495)
    label("Log hazard of re-conversion, recurrent TB -> converted (1/week)")  # Table S2 lambda32 = 0.0162 /week (RSE 21%); = 1/61.8711 ($THETA 2 MTT32) = 0.016162495
    llambda14 <- log(0.0038448829)
    label("Log baseline hazard of dropout from active TB (1/week)")  # Table S2 lambda14 = 0.00384 /week (RSE 21%); = 1/260.086 ($THETA 3 MTT14) = 0.0038448829
    llambda24 <- log(0.0014987806)
    label("Log baseline hazard of dropout from converted or recurrent TB (1/week)")  # Table S2 lambda24/34 = 0.00150 /week (RSE 18%); = 1/667.209 ($THETA 4 MTT24) = 0.0014987806. lambda34 = lambda24 (control stream MTT34 = MTT24)
    llambda25 <- log(0.00026292014)
    label("Log hazard of death from the converted state (1/week)")  # Table S2 lambda25 = 0.000263 /week (RSE 38%); = 1/3803.44 ($THETA 6 MTT25) = 0.00026292014

    # --- Weibull death hazard from active or recurrent TB --------------
    # Supplement equation 7 / control stream $DES:
    #   WB15 = scale * shape * (scale * t)^(shape - 1)
    lscale15 <- log(0.0051984779)
    label("Log Weibull scale of the death hazard from active or recurrent TB (1/week)")  # Table S2 Scale15/35 = 0.00520 /week (RSE 18%); = 1/192.364 ($THETA 5 MTT15) = 0.0051984779
    lshape15 <- log(1.96131)
    label("Log Weibull shape of the death hazard from active or recurrent TB (unitless)")  # Table S2 Shape15/35 = 1.96 (RSE 19%); $THETA 10 = 1.96131. Shape > 1 = increasing hazard over time

    # --- Covariate coefficients (log-hazard scale) ---------------------
    # Two coefficients are SHARED between two surge parameters each, with
    # opposite signs, exactly as written in the control stream $PK block.
    e_hl2_sa12pt12 <- -0.686145
    label("Coefficient of week-2 mycobacterial-load half-life on SA12 (+) and PT12 (-)")  # Table S2 betaHL2 on SA12, PT12 = -0.686 (RSE 16%); $THETA 12 = -0.686145. Centered at 0.69443 weeks
    e_ttp_pt12sw12 <- 0.442668
    label("Coefficient of baseline MGIT time-to-positivity on PT12 (-) and SW12 (+), per week of TTP")  # Table S2 betabasTTP on PT12, SW12 = 0.443 (RSE 20%); $THETA 11 = 0.442668. Centered at 217.6667 h = 9.06944 days
    e_xdr_sa12 <- -0.622792
    label("Log-hazard coefficient for XDR-TB on SA12; HR = exp(-0.623) = 0.536, a 46% reduction")  # Table S2 betaXDR-TB on SA12 = -0.623 (RSE 29%); $THETA 13 = -0.622792. Paper Results: 46% decrease (95% CI 24%-62%)
    e_sexf_lambda23 <- -0.812999
    label("Log-hazard coefficient for female sex on recurrence; HR = exp(-0.813) = 0.443 vs male")  # Table S2 betasex on lambda23 = -0.813 (RSE 41%); $THETA 17 = -0.812999
    e_mblend_lambda23 <- 0.0371081
    label("Log-hazard coefficient of end-of-treatment mycobacterial load on recurrence, per natural-log unit")  # Table S2 betaMMBLend on lambda23 = 0.0371 (RSE 41%); $THETA 18 = 0.0371081. Applied only after week 26; centered at ln(5.5726e-05)
    e_study_lambda1424 <- 0.909188
    label("Log-hazard coefficient for the C208 study on dropout; HR = exp(0.909) = 2.48 vs C209")  # Table S2 betastudy on lambda14/24/34 = 0.910 (RSE 22%); $THETA 14 = 0.909188
    e_age_lambda1424 <- -0.0230092
    label("Log-hazard coefficient of age on dropout, per year above 33 years")  # Table S2 betaage on lambda14/24/34 = -0.0230 (RSE 38%); $THETA 15 = -0.0230092
    e_baswt_lambda1525 <- -0.0838138
    label("Log-hazard coefficient of baseline body weight on death, per kg above 55 kg")  # Table S2 betabasWT on lambda15/25/35 = -0.0838 (RSE 26%); $THETA 16 = -0.0838138

    # ------------------------------------------------------------------
    # No IIV. The Lin 2024 control stream declares BSV = ETA(1) as an
    # explicit place holder and sets $OMEGA 0 FIX, so the fitted multistate
    # model is a population hazard model with no between-subject
    # variability. No eta* parameters are added.
    #
    # Residual error placeholder. The source maximises the exact multistate
    # event likelihood ($ESTIMATION METHOD=0 LIKE) on the observed
    # categorical state and has no observation-error model ($SIGMA 0 FIX in
    # the simulation variant). This nlmixr2 translation is intended for
    # forward simulation: the five state-occupancy probabilities are exposed
    # as outputs and a tiny additive residual is attached to the conversion
    # probability -- Lin 2024's primary endpoint -- so the nlmixr2 likelihood
    # machinery accepts the model. NOT from the source; see the vignette
    # Assumptions and deviations section.
    # ------------------------------------------------------------------
    addSd <- 0.001
    label("Placeholder additive residual error on the sputum-culture-conversion probability output (unitless); not from the source")
  })

  model({
    # ------------------------------------------------------------------
    # 1. Covariate-adjusted surge parameters of the conversion hazard.
    #    Control stream $PK block, verbatim:
    #      SA = THETA(7)/10000*EXP(THETA(12)*(HL2-0.69443)+THETA(13)*XDR)
    #      PT = THETA(8)*24*7*EXP(-THETA(11)*((MTTP2-217.6667)/24/7)
    #                            -THETA(12)*(HL2-0.69443))
    #      SW = THETA(9)*24*7*EXP( THETA(11)*((MTTP2-217.6667)/24/7))
    #    Note e_hl2_sa12pt12 enters SA12 with a PLUS and PT12 with a MINUS,
    #    and e_ttp_pt12sw12 enters PT12 with a MINUS and SW12 with a PLUS.
    #    The source's MTTP is in hours; TTP_MGIT_BASE is in days, so
    #    (MTTP - 217.6667)/24/7 is written as (TTP_MGIT_BASE - 9.06944)/7.
    # ------------------------------------------------------------------
    hl2Dev <- MBL_HL_WK2 - 0.69443
    ttpDev <- (TTP_MGIT_BASE - 9.06944) / 7

    sa12 <- exp(lsa12 + e_hl2_sa12pt12 * hl2Dev + e_xdr_sa12 * DIS_TB_XDR_STRICT)
    pt12 <- exp(lpt12 - e_ttp_pt12sw12 * ttpDev - e_hl2_sa12pt12 * hl2Dev)
    sw12 <- exp(lsw12 + e_ttp_pt12sw12 * ttpDev)

    # ------------------------------------------------------------------
    # 2. Transition hazards. Control stream $DES block.
    #    - lambda12: symmetric surge function (supplement equation 8).
    #    - lambda15 = lambda35: Weibull hazard (supplement equation 7). The
    #      1e-16 offset reproduces the control stream's DEL guard, which
    #      keeps (scale*t)^(shape-1) finite at t = 0 when shape < 1; it is
    #      carried over verbatim even though the estimated shape of 1.96
    #      makes the exponent positive.
    #    - The MMBLend effect on recurrence is gated on t > 26 weeks,
    #      reproducing IF(WEEK.GT.26) FLAG_MMBL = 1.
    #    - lambda14 and lambda24 have DIFFERENT baseline hazards but share a
    #      single covariate exponent; lambda34 = lambda24 and
    #      lambda35 = lambda15 by construction in the source.
    # ------------------------------------------------------------------
    dropoutCov <- e_study_lambda1424 * STUDY_C208 + e_age_lambda1424 * (AGE - 33)
    deathCov   <- e_baswt_lambda1525 * (WT - 55)

    lambda12 <- sa12 / (((t - pt12) / sw12)^2 + 1)
    lambda23 <- exp(llambda23 + e_sexf_lambda23 * SEXF +
                      e_mblend_lambda23 * (t > 26) * (log(MBL_END) - log(5.5726e-05)))
    lambda32 <- exp(llambda32)
    lambda14 <- exp(llambda14 + dropoutCov)
    lambda24 <- exp(llambda24 + dropoutCov)
    lambda34 <- lambda24
    lambda15 <- exp(lscale15) * exp(lshape15) *
      (exp(lscale15) * (t + 1e-16))^(exp(lshape15) - 1) * exp(deathCov)
    lambda25 <- exp(llambda25 + deathCov)
    lambda35 <- lambda15

    # ------------------------------------------------------------------
    # 3. Multistate ODE system (supplement equations 1-5, matching the
    #    control stream DADT(1)-DADT(5)). All patients start in active TB.
    #    s_dropout and s_death are absorbing; s_activeTb is never re-entered.
    # ------------------------------------------------------------------
    d/dt(s_activeTb) <- -s_activeTb * (lambda12 + lambda14 + lambda15)        # eq. 1
    s_activeTb(0)    <- 1

    d/dt(s_converted) <- s_activeTb * lambda12 + s_recurrentTb * lambda32 -
      s_converted * (lambda23 + lambda24 + lambda25)                         # eq. 2
    s_converted(0)    <- 0

    d/dt(s_recurrentTb) <- s_converted * lambda23 -
      s_recurrentTb * (lambda32 + lambda34 + lambda35)                       # eq. 3
    s_recurrentTb(0)    <- 0

    d/dt(s_dropout) <- s_activeTb * lambda14 + s_converted * lambda24 +
      s_recurrentTb * lambda34                                               # eq. 4
    s_dropout(0)    <- 0

    d/dt(s_death) <- s_activeTb * lambda15 + s_converted * lambda25 +
      s_recurrentTb * lambda35                                               # eq. 5
    s_death(0)    <- 0

    # ------------------------------------------------------------------
    # 4. Outputs. The five marginal state-occupancy probabilities are what
    #    Lin 2024 Figure 4 (VPC) and Figure 6 (scenario simulations) plot.
    #    haz_death reproduces the control stream $ERROR line for DV = 5,
    #    where an observed death is an exact date and so contributes the
    #    instantaneous death rate rather than a state probability.
    #    prob_state_sum is a mass-balance diagnostic: it must stay at 1.
    #    prob_scc -- sputum culture conversion, Lin 2024's primary endpoint
    #    and the one its Brier-score analysis forecasts -- carries the
    #    placeholder residual and is therefore the model's observation
    #    variable.
    # ------------------------------------------------------------------
    prob_active_tb    <- s_activeTb
    prob_scc          <- s_converted
    prob_recurrent_tb <- s_recurrentTb
    prob_dropout      <- s_dropout
    prob_death        <- s_death

    haz_death <- s_activeTb * lambda15 + s_converted * lambda25 +
      s_recurrentTb * lambda35
    prob_state_sum <- s_activeTb + s_converted + s_recurrentTb + s_dropout + s_death

    prob_scc ~ add(addSd)
  })
}
