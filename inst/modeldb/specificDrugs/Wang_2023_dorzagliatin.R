Wang_2023_dorzagliatin <- function() {
  description <- paste(
    "Two-compartment population PK model for oral dorzagliatin (HMS5552,",
    "a first-in-class dual-acting glucokinase activator approved in China",
    "for type 2 diabetes mellitus) in healthy subjects and patients with",
    "T2DM (Wang 2023; N = 1062 subjects, 7686 concentrations pooled from",
    "six clinical trials). Disposition is a two-compartment model with",
    "sequential zero-order then first-order oral absorption (dose enters",
    "the depot over a zero-order window of duration D1 = 0.418 h, then is",
    "absorbed first-order at ka = 3.29 h^-1) and first-order elimination",
    "from the central compartment. Typical values for a 69 kg, 55-year-old",
    "male with AST 18 U/L enrolled in the later-phase studies are CL/F =",
    "10.4 L/h, Vc/F = 80.6 L, Q/F = 3.02 L/h and Vp/F = 26.5 L. CL/F",
    "carries three power covariate effects plus one study indicator:",
    "(WT/69)^0.255, (AST/18)^-0.103, (AGE/55)^-0.135 and a multiplicative",
    "exp(0.203) = 1.23 uplift for the three early-phase studies",
    "(HMM0102 / HMM0103 / HMM0110), which the paper attributes to a",
    "formulation / bioavailability difference that could not be resolved",
    "as an F effect. Vc/F carries (WT/69)^0.553 and a multiplicative",
    "exp(-0.170) = 0.843 factor in females (15.7% lower than males). D1",
    "carries exp(0.816) = 2.26 when food is consumed at least 1 h after",
    "dosing rather than 0.5 h after dosing. Because sex acts only on Vc/F",
    "and the meal-timing effect only on D1, neither changes steady-state",
    "AUCtau (the paper reports 0.01% and 0% respectively) while both shift",
    "Cmax,ss and Cmin,ss. Inter-individual variability is a correlated",
    "CL/F-Vc/F block (approximate CV 22.5% and 14.9%; covariance 0.0181)",
    "plus independent terms on Vp/F (48.8%) and D1 (82.8%); ka and Q/F",
    "carry no IIV. Residual error is combined additive plus proportional",
    "(109 ng/mL and 32.9%).",
    sep = " "
  )
  reference <- paste(
    "Wang K, Feng L, Zhang J, Zou Q, Xu F, Sun Z, Tang F, Chen L. (2023).",
    "Population pharmacokinetic analysis of dorzagliatin in healthy",
    "subjects and patients with type 2 diabetes mellitus.",
    "Clin Pharmacokinet 62(10):1419-1430.",
    "doi:10.1007/s40262-023-01286-8",
    sep = " "
  )
  vignette <- "Wang_2023_dorzagliatin"
  units    <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Baseline, time-fixed. Overall population median 69.0 kg",
        "(range 40.0-110; Wang 2023 Table 2 'TBW, kg' Overall column),",
        "which is the normalization constant used by both covariate",
        "relationships. Retained on CL/F and Vc/F in the final model as",
        "power effects (WT/69)^0.255 and (WT/69)^0.553 -- the paper prints",
        "these in the algebraically identical exp(theta * log(WT/69))",
        "form (Wang 2023 Table 3 footnote equations). Note that neither",
        "exponent is the allometric 0.75 / 1.0 pair; both were estimated",
        "(RSE 19.4% and 9.94%). Wang 2023 Discussion reports the",
        "resulting 10th-90th percentile spread (55 kg and 83 kg) as",
        "-5.63% to +4.83% on CL/F and -11.8% to +10.8% on Vc/F."
      ),
      source_name        = "TBW"
    ),
    AGE = list(
      description        = "Age at baseline",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Baseline, time-fixed. Overall population median 55.0 years",
        "(range 19.0-74.2; Wang 2023 Table 2 'AGE, year' Overall column),",
        "used as the normalization constant. Retained on CL/F as the",
        "power effect (AGE/55)^-0.135, i.e. older subjects have lower",
        "CL/F and therefore higher exposure. Wang 2023 Discussion",
        "reports +4.06% (41 years) to -2.44% (66 years) on CL/F across",
        "the 10th-90th percentiles. No subject was younger than 19 or",
        "older than 74.2 years, so the relationship is not established",
        "in pediatric or very elderly subjects."
      ),
      source_name        = "AGE"
    ),
    AST = list(
      description        = "Baseline aspartate aminotransferase",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Baseline, time-fixed. Overall population median 18.0 U/L",
        "(range 8.00-74.0; Wang 2023 Table 2 'AST, U/L' Overall column),",
        "used as the normalization constant. Reported by the paper in",
        "U/L, matching the canonical SI unit (the source uses U/L and",
        "IU/L interchangeably; the supplement quotes the AST upper limit",
        "of normal as 40 IU/L). Retained on CL/F as the power effect",
        "(AST/18)^-0.103: higher AST (poorer hepatocellular integrity)",
        "gives lower CL/F and higher exposure. Wang 2023 Discussion",
        "reports +3.41% (13 U/L) to -4.45% (28 U/L) on CL/F across the",
        "10th-90th percentiles. AST is the only retained hepatic marker;",
        "ALT and total bilirubin were screened and dropped."
      ),
      source_name        = "AST"
    ),
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "Time-fixed per subject. The source column GEND is coded 1 =",
        "female, 0 = male (Wang 2023 Table 3 footnote: 'GEND: 1 for",
        "female, 0 for male'), which is exactly the canonical SEXF",
        "orientation -- no value transformation is required. Cohort",
        "397 / 1062 female (37.4%) per Wang 2023 Table 2. Retained on",
        "Vc/F only as a multiplicative exp(-0.170) = 0.843 factor,",
        "i.e. females have 15.7% lower Vc/F than males (Wang 2023",
        "Discussion). Because sex does not act on CL/F, steady-state",
        "AUCtau is essentially unchanged between sexes (0.01% per",
        "Wang 2023 Section 3.3) while Cmax,ss is 9.09% higher and",
        "Cmin,ss 11.6% lower in females."
      ),
      source_name        = "GEND"
    ),
    STUDY_DORZA_EARLY = list(
      description        = "Early-phase dorzagliatin study indicator (HMM0102 / HMM0103 / HMM0110)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (later-phase studies HMM0201 / HMM0301 / HMM0302)",
      notes              = paste(
        "Time-fixed per subject; set from the trial identifier. 1 = the",
        "subject was enrolled in study HMM0102 (NCT02077452), HMM0103",
        "(NCT02386982) or HMM0110 (NCT04324424); 0 = HMM0201",
        "(NCT02561338), HMM0301 (NCT03173391) or HMM0302 (NCT03141073).",
        "Wang 2023 Table 3 footnote: 'Study: 1 for study 102/103/110, 0",
        "for others'. Enters CL/F as exp(0.203 * STUDY_DORZA_EARLY) =",
        "a 1.23-fold (22.5% higher) apparent clearance in the",
        "early-phase studies. Wang 2023 Discussion paragraph 2 explains",
        "the mechanism: those studies showed systematically lower",
        "exposure, plausibly a formulation difference (studies 301 and",
        "302 used near-commercial preparations), but 'the model was",
        "unable to correct for prediction bias by bioavailability' in",
        "the sparsely-sampled later studies, so the between-trial",
        "difference was absorbed into CL/F instead of F."
      ),
      source_name        = "Study"
    ),
    MEAL_DELAY_1H = list(
      description        = "Delayed-meal indicator: food consumed at least 1 h after dosing rather than 0.5 h after dosing",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (food consumed 0.5 h after dosing)",
      notes              = paste(
        "Per-dose-record indicator. 1 = the subject started eating at",
        "least 1 h after taking the dose; 0 = the subject started eating",
        "0.5 h after taking the dose (Wang 2023 Table 3 footnote:",
        "'FOOD: 0 and 1 for the time of food consumption after 0.5 h and",
        ">= 1 h of drug administration, respectively'). Both levels are",
        "fed states -- this covariate encodes the dose-to-meal interval,",
        "not the presence or the fat content of a meal, so it is distinct",
        "from FED, FED_HIGHFAT and FED_LOWFAT. Retained on D1 only, as",
        "the multiplicative factor exp(0.816) = 2.26: waiting at least",
        "1 h before eating more than doubles the zero-order absorption",
        "window (0.418 h -> 0.945 h). Because D1 does not affect CL/F,",
        "steady-state AUCtau is unchanged and Cmax,ss / Cmin,ss shift by",
        "only -1.26% / +2.71% (Wang 2023 Section 3.3). Note that Wang",
        "2023 Table 3 mislabels this row as 'D1,FOOD -- Influence of",
        "food consumption time on CL/F'; the printed equation and the",
        "parameter name both place the effect on D1, and the Discussion",
        "confirms 'patients with time to food consumption of >= 1 h",
        "after drug administration had higher D1'."
      ),
      source_name        = "FOOD"
    )
  )

  covariatesDataExcluded <- list(
    BMI = list(
      description        = "Baseline body mass index",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened in the stepwise forward-inclusion / backward-",
        "elimination covariate search (Wang 2023 Section 2.4) but not",
        "retained in the final model. Overall median 25.2 kg/m^2",
        "(range 18.2-34.9) per Table 2. Wang 2023 Section 3.4 and",
        "Table S1 report the geometric-mean exposure difference across",
        "BMI quartiles as only -10% to +5.57% relative to the overall",
        "population, driven by correlated covariates rather than a",
        "retained BMI effect."
      ),
      source_name        = "BMI"
    ),
    BSA = list(
      description        = "Baseline body surface area",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened (Wang 2023 Section 2.4 covariate list) but not",
        "retained; body weight was the size descriptor kept on both",
        "CL/F and Vc/F. Not tabulated in Wang 2023 Table 2."
      ),
      source_name        = "BSA"
    ),
    ALB = list(
      description        = "Baseline serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened but not retained. Overall median 46.2 g/L (range",
        "33.2-57.0) per Wang 2023 Table 2; already in canonical SI",
        "units, so no conversion is needed."
      ),
      source_name        = "ALB"
    ),
    ALT = list(
      description        = "Baseline alanine aminotransferase",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened but not retained; AST was the hepatic marker kept on",
        "CL/F. Overall median 19.0 U/L (range 2.50-110) per Wang 2023",
        "Table 2. Wang 2023 Section 3.2 states explicitly that ALT",
        "'had no significant effect on the pharmacokinetics of",
        "dorzagliatin'."
      ),
      source_name        = "ALT"
    ),
    CREAT = list(
      description        = "Baseline serum creatinine",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened but not retained. Overall median 67.0 umol/L (range",
        "35.0-943) per Wang 2023 Table 2; the extreme upper end comes",
        "from study HMM0110, the dedicated renal-impairment study",
        "(median 235 umol/L in that cohort). Already in canonical SI",
        "units."
      ),
      source_name        = "CR"
    ),
    CRCL = list(
      description        = "Baseline creatinine clearance",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened but not retained. Not tabulated in Wang 2023 Table 2",
        "(the reported renal marker is serum creatinine). Wang 2023",
        "Section 3.2 states creatinine clearance 'had no significant",
        "effect'; the Discussion adds that a two-part renal-impairment",
        "study found only limited exposure differences between ESRD",
        "subjects and healthy volunteers (geometric mean ratio 0.81 for",
        "Cmax and 1.11 for AUC), so no dose adjustment is needed at any",
        "stage of renal impairment. The paper does not state whether",
        "the screened column was BSA-normalized."
      ),
      source_name        = "CRCL"
    ),
    TBILI = list(
      description        = "Baseline total bilirubin",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened but not retained. Overall median 10.0 umol/L (range",
        "1.25-40.0) per Wang 2023 Table 2; already in canonical SI",
        "units. Used together with AST to assign the FDA hepatic-",
        "function categories in Wang 2023 Table S3 (the supplement",
        "quotes upper limits of normal of 20.5 umol/L for total",
        "bilirubin and 40 IU/L for AST), but it carries no retained",
        "effect on any PK parameter."
      ),
      source_name        = "TBIL"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 1062L,
    n_studies       = 6L,
    n_observations  = 7686L,
    age_range       = "19.0-74.2 years",
    age_median      = "55.0 years",
    weight_range    = "40.0-110 kg",
    weight_median   = "69.0 kg",
    bmi_median      = "25.2 kg/m^2 (range 18.2-34.9)",
    sex_female_pct  = 37.4,
    race_ethnicity  = c(Asian = 100),
    disease_state   = paste(
      "Pooled cohort of healthy volunteers and adult patients with type 2",
      "diabetes mellitus from six trials (Wang 2023 Table 1): two phase I",
      "studies (HMM0102 / NCT02077452 multiple ascending dose;",
      "HMM0103 / NCT02386982 28-day beta-cell-function study), one phase",
      "II dose-ranging study (HMM0201 / NCT02561338), two phase III",
      "studies (HMM0301 / NCT03173391 monotherapy in drug-naive patients;",
      "HMM0302 / NCT03141073 add-on to metformin), and one open-label",
      "renal-impairment study with matched healthy volunteers",
      "(HMM0110 / NCT04324424)."
    ),
    dose_range      = paste(
      "25-400 mg/day orally: 25, 50, 100, 150 and 200 mg twice daily",
      "(HMM0102); 75 mg once or twice daily and 100 mg once daily",
      "(HMM0103, HMM0201, HMM0301, HMM0302); 25 mg twice daily",
      "(HMM0110). The approved and pivotal-trial regimen is 75 mg",
      "twice daily."
    ),
    regions         = "China (all six trials were conducted by Hua Medicine in Chinese subjects)",
    hepatic_function = paste(
      "Per the FDA hepatic-impairment categories derived from AST and",
      "total bilirubin (Wang 2023 Table S3): normal 974 / 92.2%, mild",
      "77 / 7.29%, moderate 5 / 0.473%, severe 0. Baseline AST median",
      "18.0 U/L (8.00-74.0), ALT 19.0 U/L (2.50-110), albumin 46.2 g/L",
      "(33.2-57.0), total bilirubin 10.0 umol/L (1.25-40.0)."
    ),
    renal_function  = paste(
      "Baseline serum creatinine median 67.0 umol/L (35.0-943); the wide",
      "upper range comes from the dedicated renal-impairment study",
      "HMM0110 (cohort median 235 umol/L), which enrolled subjects",
      "through end-stage renal disease. Wang 2023 Table S4 groups the",
      "analysis population as normal 936 and mild 109 by eGFR."
    ),
    notes           = paste(
      "Wang 2023 Table 2 summarises baseline demographics by study and",
      "overall; Table 1 lists the study designs, dose regimens and PK",
      "sampling schedules. Parameter estimates implemented here are the",
      "final-model 'Estimate (%RSE)' column of Wang 2023 Table 3, fit in",
      "NONMEM 7.4.3 with FOCE-with-interaction and confirmed by a",
      "1000-replicate nonparametric bootstrap whose medians agree with",
      "the point estimates to within rounding. Race is not tabulated;",
      "all six trials were conducted in China, so the cohort is recorded",
      "here as Asian."
    )
  )

  ini({
    # -----------------------------------------------------------------------
    # Structural PK -- Wang 2023 Table 3, 'Final model Estimate (%RSE)'
    # column. All disposition parameters are apparent (i.e. /F): the paper
    # analysed oral data only and did not estimate an absolute
    # bioavailability, so F is implicitly 1 and no lfdepot term is used.
    # Typical values correspond to a 69 kg, 55-year-old male with AST
    # 18 U/L enrolled in the later-phase studies (Wang 2023 Section 3.2).
    # Units are hour / L / L/h; concentration is converted to ng/mL inside
    # model() (see the observation comment there).
    # -----------------------------------------------------------------------
    lcl <- log(10.4)  ; label("Apparent total clearance CL/F (L/h)")                        # Wang 2023 Table 3 CL/F = 10.4 L/h (RSE 0.923%, 95% CI 10.2-10.6; bootstrap median 10.4)
    lvc <- log(80.6)  ; label("Apparent central volume of distribution Vc/F (L)")           # Wang 2023 Table 3 Vc/F = 80.6 L (RSE 1.23%, 95% CI 78.7-82.6; bootstrap median 80.6)
    lq  <- log(3.02)  ; label("Apparent inter-compartmental clearance Q/F (L/h)")           # Wang 2023 Table 3 Q/F = 3.02 L/h (RSE 9.70%, 95% CI 2.50-3.66; bootstrap median 3.01)
    lvp <- log(26.5)  ; label("Apparent peripheral volume of distribution Vp/F (L)")        # Wang 2023 Table 3 Vp/F = 26.5 L (RSE 6.95%, 95% CI 23.1-30.4; bootstrap median 26.7)
    lka <- log(3.29)  ; label("First-order absorption rate constant ka (1/h)")              # Wang 2023 Table 3 Ka = 3.29 1/h (RSE 5.30%, 95% CI 2.96-3.65; bootstrap median 3.29)
    ld1 <- log(0.418) ; label("Zero-order absorption input duration D1 (h)")                # Wang 2023 Table 3 D1 = 0.418 h (RSE 7.49%, 95% CI 0.361-0.484; bootstrap median 0.415)

    # -----------------------------------------------------------------------
    # Covariate effects. Wang 2023 prints the three final-model covariate
    # relations immediately below Table 3 as:
    #
    #   CL/F = 10.4 * exp(0.203 * Study + 0.255 * log(BW/69)
    #                     - 0.103 * log(AST/18) - 0.135 * log(AGE/55) + eta_CL)
    #   Vc/F = 80.6 * exp(0.553 * log(BW/69) - 0.170 * GEND + eta_Vc)
    #   D1   = 0.418 * exp(0.816 * FOOD + eta_D1)
    #
    # The continuous terms exp(theta * log(cov/ref)) are written below in the
    # algebraically identical power form (cov/ref)^theta, which is the
    # nlmixr2lib convention. The binary terms stay in exp() form.
    #
    # Table 3 reports the three binary/indicator effects already
    # exponentiated, which is how the two representations reconcile:
    #   CL_STUDY  = 1.23  = exp(0.203)   -> +22.5% CL/F in the early studies
    #   Vc,GEND   = 0.843 = exp(-0.170)  -> -15.7% Vc/F in females
    #   D1,FOOD   = 2.26  = exp(0.816)   -> +126% D1 when eating >= 1 h post-dose
    # while the continuous power exponents are tabulated on the log scale
    # directly (CL_BW 0.255, Vc,BW 0.553, CL_AST -0.103, CL_AGE -0.135).
    # Every one of these coefficients was checked by reproducing the paper's
    # own reported percentage changes -- see the vignette source-trace table.
    # -----------------------------------------------------------------------
    e_study_dorza_early_cl <- 0.203  ; label("Early-phase-study effect on CL/F (log-scale indicator coefficient)")           # Wang 2023 Table 3 footnote equation '0.203 x Study'; Table 3 CL_STUDY = 1.23 = exp(0.203) (RSE 2.39%, 95% CI 1.17-1.28)
    e_wt_cl                <- 0.255  ; label("Power effect of body weight on CL/F (unitless, reference 69 kg)")             # Wang 2023 Table 3 CL_BW = 0.255 (RSE 19.4%, 95% CI 0.158-0.352); footnote equation '0.255 x log(BW/69)'
    e_ast_cl               <- -0.103 ; label("Power effect of AST on CL/F (unitless, reference 18 U/L)")                    # Wang 2023 Table 3 CL_AST = -0.103 (RSE 21.3%, 95% CI -0.146 to -0.0601); footnote equation '-0.103 x log(AST/18)'
    e_age_cl               <- -0.135 ; label("Power effect of age on CL/F (unitless, reference 55 years)")                  # Wang 2023 Table 3 CL_AGE = -0.135 (RSE 28.5%, 95% CI -0.211 to -0.0599); footnote equation '-0.135 x log(AGE/55)'
    e_wt_vc                <- 0.553  ; label("Power effect of body weight on Vc/F (unitless, reference 69 kg)")             # Wang 2023 Table 3 Vc,BW = 0.553 (RSE 9.94%, 95% CI 0.445-0.661); footnote equation '0.553 x log(BW/69)'
    # NOTE on e_sexf_vc: the paper quotes this effect three ways that do not
    # round to one another -- equation coefficient -0.170, Table 3
    # Vc,GEND = 0.843, and Discussion "females were 15.7% lower". But
    # exp(-0.170) = 0.84366 (rounds to 0.844, not 0.843) and
    # log(0.843) = -0.170788. Both are roundings of the same underlying
    # estimate (~-0.1708) at different precisions. The equation coefficient is
    # used here per the printed-equation-governs convention; the resulting
    # female/male Vc/F ratio is 0.8437 (-15.63%) rather than the quoted -15.7%.
    # See the vignette Errata for the measured consequence on Cmax,ss/Cmin,ss
    # (all differences < 0.07 percentage points).
    e_sexf_vc              <- -0.170 ; label("Female-sex effect on Vc/F (log-scale indicator coefficient)")                 # Wang 2023 Table 3 footnote equation '-0.170 x GEND'; cf. Table 3 Vc,GEND = 0.843 (RSE 1.74%, 95% CI 0.815-0.873)
    e_meal_delay_1h_d1     <- 0.816  ; label("Delayed-meal effect on D1 (log-scale indicator coefficient)")                 # Wang 2023 Table 3 footnote equation '0.816 x FOOD'; Table 3 D1,FOOD = 2.26 = exp(0.816) (RSE 12.8%, 95% CI 1.76-2.91)

    # -----------------------------------------------------------------------
    # Inter-individual variability. Wang 2023 Table 3 reports IIV_CL 22.5%,
    # IIV_Vc 14.9%, IIV_Vp 48.8% and IIV_D1 82.8% with the footnote "IIV for
    # CL, Vc, Vp, D1, and proportional error are reported as approximate
    # CV%". "Approximate" identifies the sqrt-of-variance convention
    # CV% = 100 * sqrt(omega^2), i.e. omega = CV/100, rather than the exact
    # log-normal CV% = 100 * sqrt(exp(omega^2) - 1); the exact form would not
    # be described as an approximation. The two readings differ negligibly
    # here in any case (omega_CL 0.2250 vs 0.2222).
    #
    # Variances entered below are therefore (CV/100)^2:
    #   var(CL) = 0.225^2 = 0.050625      var(Vc) = 0.149^2 = 0.022201
    #   var(Vp) = 0.488^2 = 0.238144      var(D1) = 0.828^2 = 0.685584
    # The CL-Vc covariance is taken directly from Table 3 Cor_CL,Vc =
    # 0.0181, which is labelled "Covariance (CL, Vc)" and is already on the
    # omega scale. It implies a correlation of
    # 0.0181 / (0.225 * 0.149) = 0.540, and satisfies the positive-definite
    # bound |cov| < sqrt(var_CL * var_Vc) = 0.0335.
    #
    # ka and Q/F have no IIV term in Wang 2023 Table 3 and none is invented
    # here (the model is encoded exactly as published).
    # -----------------------------------------------------------------------
    etalcl + etalvc ~ c(0.050625,
                        0.0181, 0.022201)                                                                                   # Wang 2023 Table 3: IIV_CL 22.5% CV -> var 0.050625; Cor_CL,Vc = 0.0181 (RSE 10.9%); IIV_Vc 14.9% CV -> var 0.022201
    etalvp ~ 0.238144                                                                                                       # Wang 2023 Table 3 IIV_Vp = 48.8% CV (RSE 9.25%, 95% CI 38.9-57.0) -> var 0.488^2
    etald1 ~ 0.685584                                                                                                       # Wang 2023 Table 3 IIV_D1 = 82.8% CV (RSE 8.16%, 95% CI 68.3-95.2) -> var 0.828^2

    # -----------------------------------------------------------------------
    # Residual error: combined additive plus proportional (Wang 2023 Section
    # 2.3, "The final residual error model included additive and
    # proportional error terms"). Both Table 3 rows are reported on the SD
    # scale in their own units -- 'sigma1 (%) Proportional error 32.9' and
    # 'sigma2 (ng/mL) Additive error 109' -- so they are used directly as
    # propSd and addSd rather than being square-rooted.
    # -----------------------------------------------------------------------
    propSd <- 0.329 ; label("Proportional residual error (fraction)")     # Wang 2023 Table 3 sigma1 = 32.9% (RSE 1.34%, 95% CI 32.0-33.7; bootstrap median 32.8)
    addSd  <- 109   ; label("Additive residual error (ng/mL)")            # Wang 2023 Table 3 sigma2 = 109 ng/mL (RSE 1.38%, 95% CI 106-112; bootstrap median 109)
  })

  model({
    # ---------------------------------------------------------------------
    # 1. Individual PK parameters.
    #
    # CL/F carries the early-study indicator plus three power covariate
    # effects; Vc/F carries body weight and female sex; Vp/F, Q/F and ka
    # are covariate-free. D1 carries the meal-timing indicator.
    #
    # Reference (normalization) values are the overall population medians
    # from Wang 2023 Table 2: 69 kg, 18 U/L AST, 55 years. At the reference
    # covariate vector with STUDY_DORZA_EARLY = 0, SEXF = 0 and
    # MEAL_DELAY_1H = 0 every multiplier collapses to 1 and the parameters
    # reduce to the Table 3 typical values.
    # ---------------------------------------------------------------------
    cl <- exp(lcl + etalcl) *
      exp(e_study_dorza_early_cl * STUDY_DORZA_EARLY) *
      (WT / 69)^e_wt_cl *
      (AST / 18)^e_ast_cl *
      (AGE / 55)^e_age_cl
    vc <- exp(lvc + etalvc) *
      exp(e_sexf_vc * SEXF) *
      (WT / 69)^e_wt_vc
    vp <- exp(lvp + etalvp)
    q  <- exp(lq)
    ka <- exp(lka)
    d1 <- exp(ld1 + etald1) * exp(e_meal_delay_1h_d1 * MEAL_DELAY_1H)

    # 2. Micro-constants for the two-compartment system.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ---------------------------------------------------------------------
    # 3. Two-compartment ODE system with sequential zero-order then
    #    first-order absorption (Wang 2023 Section 2.3 / 3.2). The dose is
    #    delivered into `depot` at a constant rate over the zero-order
    #    window of duration D1, and `depot` then drains into `central`
    #    first-order at ka. Dose records must therefore carry rate = -2 so
    #    that rxode2 uses the model's dur(depot); a plain bolus would
    #    collapse the zero-order phase and bias Cmax.
    # ---------------------------------------------------------------------
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    dur(depot) <- d1

    # ---------------------------------------------------------------------
    # 4. Observation. Dose is in mg and vc in L, so central / vc is mg/L
    #    (= ug/mL); multiplying by 1000 gives ng/mL, the assay unit used
    #    throughout Wang 2023 (Table 3 additive residual error 109 ng/mL;
    #    Table S1-S4 steady-state exposures ~1088 ng/mL Cmax,ss and
    #    ~7153 ng*h/mL AUCtau,ss at 75 mg twice daily).
    # ---------------------------------------------------------------------
    Cc <- central / vc * 1000
    Cc ~ add(addSd) + prop(propSd)
  })
}
