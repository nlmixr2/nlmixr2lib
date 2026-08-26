Gosselin_2015_motesanib <- function() {
  description <- "Joint parent plus metabolite population PK model for motesanib and its active metabolite M4 (indoline lactam) in patients with advanced solid tumors, pooled from 8 phase 1 to phase 3 oral-dosing studies. Motesanib disposition is two-compartment with first-order absorption, an absorption lag time, and a high-fat-meal food effect on both the absorption rate constant and the lag time; apparent clearance carries a near-linear power effect of serum albumin plus a multiplicative sex effect, and apparent central volume carries power effects of albumin, alkaline phosphatase and body weight. M4 is described by a one-compartment model whose formation flux is the entire apparent elimination flux of motesanib, so its apparent clearance and volume are scaled by FM4 (the product of the oral absorbed fraction, the fraction biotransformed to M4, and the parent-to-metabolite molecular-weight ratio); M4 apparent clearance carries Asian-race and once-daily-dosing effects and M4 apparent volume carries sex and once-daily-dosing effects. Inter-occasion variability over two occasions (the first week of treatment versus the remaining weeks) is carried on motesanib apparent clearance and on both M4 disposition parameters."
  reference <- paste(
    "Gosselin N. H., Mouksassi M.-S., Lu J.-F., Hsu C.-P. (2015).",
    "Population pharmacokinetic modeling of motesanib and its active",
    "metabolite, M4, in cancer patients.",
    "Clinical Pharmacology in Drug Development 4(6):463-472.",
    "doi:10.1002/cpdd.196.",
    sep = " "
  )
  vignette <- "Gosselin_2015_motesanib"

  # Gosselin 2015 reports plasma concentrations in ug/L (Table 4: motesanib
  # Cmax 493.4 ug/L after 125 mg once daily; 125 mg / 240 L = 0.52 mg/L =
  # 520 ug/L) and the fixed additive residual errors in ug/L. Clearances are
  # L/h and volumes are L, so an amount in mg over a volume in L gives mg/L;
  # the model therefore carries an explicit ugPerMg factor on the two
  # observation lines so that Cc and Cc_m4 come out in the paper's ug/L.
  units <- list(time = "h", dosing = "mg", concentration = "ug/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Amounts are APPARENT (divided by the relevant
  # bioavailability term) because every disposition parameter in the paper is
  # an apparent value: CL/F, Vc/F, Q/F, Vp/F for motesanib and CLM4/FM4,
  # VM4/FM4 for M4.
  compartmentData <- list(
    depot = list(
      analyte = "motesanib", units = "mg (apparent, i.e. amount/F)",
      specimen = "administration site", verified = FALSE
    ),
    central = list(
      analyte = "motesanib", units = "mg (apparent, i.e. amount/F)",
      specimen = "plasma", verified = FALSE
    ),
    peripheral1 = list(
      analyte = "motesanib", units = "mg (apparent, i.e. amount/F)",
      specimen = "plasma", verified = FALSE
    ),
    central_m4 = list(
      analyte = "M4", units = "mg motesanib-equivalents (apparent, i.e. amount/FM4)",
      specimen = "plasma", verified = FALSE
    )
  )

  covariateData <- list(
    ALB = list(
      description        = "Baseline serum total albumin, the most influential covariate on motesanib CL/F and also a covariate on Vc/F",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Entered as a power function centred at 39 g/L on both CL/F and Vc/F",
        "(Gosselin 2015 Table 3 'ALB on CL/F = (ALB/39)^0.971' and 'ALB on Vc/F",
        "= (ALB/39)^1.66'). Methods 'Population PK Modeling of Motesanib':",
        "'Continuous covariates were included in the structural model with a",
        "power function centering at median values'. The centring value 39 g/L",
        "equals the cohort median in Table 2 (39.0 g/L, 5th-95th percentiles",
        "30.0-46.0) and is restated in the Results text ('These later values",
        "were based on central levels for ALB and ALP (ie, 39 g/L and 101 U/L,",
        "respectively)'). Observed albumin range was 22-48 g/L in men and",
        "24-55 g/L in women (Results, Figure 2 legend); missing albumin for 3",
        "patients was imputed to the median before modelling (Methods 'Data",
        "Assembly'). The paper reports the exponent 0.971 on CL/F as 'close to",
        "1', i.e. the CL/F-albumin relationship is nearly proportional, and",
        "notes that predicted clearance is at least 22% lower in patients with",
        "hypoalbuminemia (< 34 g/L) than in patients with albumin 34-54 g/L."
      ),
      source_name        = "ALB"
    ),
    ALP = list(
      description        = "Baseline serum alkaline phosphatase, a covariate on motesanib Vc/F",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Entered as a power function centred at 101 U/L on Vc/F (Gosselin 2015",
        "Table 3 'ALP on Vc/F = (ALP/101)^-0.217'). The negative exponent",
        "matches the exploratory finding of 'a negative trend between EBEs of",
        "Vc/F and ALP' (Results 'Population PK Modeling of Motesanib') and the",
        "bootstrap median printed alongside it (-0.218). The centring value",
        "101 U/L is stated verbatim in the Results text; it differs from the",
        "Table 2 pooled median of 106.0 U/L because Table 2 describes the",
        "combined 451-patient motesanib + M4 cohort while the motesanib model",
        "was fitted to 445 patients. ALP was missing for 1 patient and was",
        "imputed to the median before modelling (Methods 'Data Assembly')."
      ),
      source_name        = "ALP"
    ),
    WT = list(
      description        = "Baseline body weight, a covariate on motesanib Vc/F",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Entered as a power function centred at 68.3 kg on Vc/F (Gosselin 2015",
        "Table 3 'WT on Vc/F = (WT/68.3)^0.612'). This is an ESTIMATED",
        "exponent, not a fixed allometric 1: the paper's stepwise covariate",
        "search retained it at 0.799 in step 2 and it settled at 0.612 in the",
        "final joint fit (Supplement Table 1 'STEP 2 +WT on Vc/F, estimate",
        "0.799'). The centring value 68.3 kg is close to but not identical",
        "with the Table 2 pooled median of 67.41 kg (5th-95th percentiles",
        "45.0-103.8) because Table 2 describes the combined 451-patient cohort",
        "while the motesanib model was fitted to 445 patients. Weight was",
        "tested on CL/F as well but was NOT retained (Results: the formal",
        "covariate analysis tested 'CL/F: sex, race (Asian versus non-Asian),",
        "weight, age, ALP, AST, albumin (ALB), serum creatinine, and total",
        "bilirubin')."
      ),
      source_name        = "WT"
    ),
    SEXF = list(
      description        = "Female sex indicator; a covariate on motesanib CL/F and on M4 VM4/FM4",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "1 = female, 0 = male. Applied as a multiplicative factor on the male",
        "typical value in both places: motesanib CL/F is 0.882-fold in women",
        "(Gosselin 2015 Table 3 'Sex on CL/F = 0.882 for women'; the male",
        "typical value 62.7 L/h times 0.882 gives 55.3 L/h, exactly the pair",
        "'62.7 vs 55.3 L/h' quoted in the Discussion), and M4 VM4/FM4 is",
        "0.57-fold in women (Table 3 'Sex on VM4/FM4 = 0.57 for women'; the",
        "Discussion restates this as 'a volume approximately 1.8-fold higher",
        "for male patients', and 1/0.57 = 1.75). Note that the bootstrap",
        "median printed for the CL/F sex factor in Table 3 reads '(0.0882)',",
        "an obvious decimal-point typographical error for 0.882: 62.7 * 0.0882",
        "would be 5.53 L/h, which contradicts the Discussion's 55.3 L/h and",
        "the reported female terminal half-life of 6.20 h. The cohort was 232",
        "men / 219 women (Table 2)."
      ),
      source_name        = "SEX"
    ),
    FED_HIGHFAT = list(
      description        = "High-fat, high-calorie meal indicator at the dosing occasion; the paper's 'diet' covariate on the motesanib absorption rate constant and lag time",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted)",
      notes              = paste(
        "1 = dose taken 5 minutes after a standardized high-fat, high-caloric",
        "breakfast; 0 = fasted. Gosselin 2015 Methods 'Clinical Studies': 'In",
        "most studies, motesanib was taken under fasted conditions (no food or",
        "liquids, except water, for 1 or 2 hours before motesanib intake). In",
        "1 study, food effect on motesanib PK was assessed in 10 patients with",
        "intensive PK samples collected after dosing, occurring 5 minutes",
        "after eating a standardized high-fat, high-caloric breakfast.' The",
        "canonical FED_HIGHFAT (rather than the more general FED) is used",
        "because every fed record in the analysis dataset came from that",
        "high-fat-breakfast substudy; the paper itself labels the covariate",
        "generically as 'diet' or 'food status (fasted vs fed status)'.",
        "Time-VARYING within a subject: the food-effect substudy dosed the",
        "same patients fed on one occasion and fasted on others. Applied as a",
        "multiplicative factor on both absorption parameters (Table 3 'Diet on",
        "lag time = x 1.079' and 'Diet on Ka = x 0.0245'), so a fed dose has a",
        "7.9% longer lag and an absorption rate constant 40-fold lower",
        "(9.84 -> 0.241 1/h). Dosing was assumed to be fasted when diet status",
        "was missing (Methods 'Data Assembly')."
      ),
      source_name        = "DIET"
    ),
    RACE_ASIAN = list(
      description        = "Asian race indicator; a covariate on M4 apparent clearance",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Asian)",
      notes              = paste(
        "1 = Asian, 0 = non-Asian. Applied as a multiplicative factor on the",
        "non-Asian typical value of CLM4/FM4 (Gosselin 2015 Table 3 'Asian on",
        "CLM4/FM4 = x 1.4 for Asian'), i.e. M4 apparent clearance is 1.4-fold",
        "higher in Asian patients; the Abstract calls this 'faster elimination",
        "of M4 in Asian patients'. The dichotomy tested was explicitly 'race",
        "(Asian versus non-Asian)' (Results). The cohort contained 142 Asians",
        "(31.5%), of whom 102 (72%) were Japanese (Results 'Patient",
        "Demographics'). Race was also tested on motesanib CL/F and Vc/F and",
        "on M4 VM4/FM4 but was NOT retained in those places."
      ),
      source_name        = "ASIAN"
    ),
    REGI_QD = list(
      description        = "Once-daily dosing-regimen indicator; a covariate on both M4 disposition parameters, with twice-daily dosing as the reference",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (twice-daily, every-12-hour dosing)",
      notes              = paste(
        "1 = the patient's motesanib regimen is once daily (every 24 h),",
        "0 = twice daily (every 12 h). Gosselin 2015 makes TWICE-daily the",
        "reference and reports both effects 'for once-daily dose' (Table 3",
        "'Dosing interval on CLM4/FM4 = x 1.52 for once-daily dose' and",
        "'Dosing interval on VM4/FM4 = x 0.45 for once-daily dose'), which is",
        "why the once-daily indicator REGI_QD is used here rather than the",
        "complementary canonical REGI_BID. The two are mutually exclusive and",
        "exhaustive in this analysis (REGI_QD = 1 - REGI_BID); the paper",
        "pooled only q24h and q12h arms. This reference orientation is",
        "confirmed arithmetically by the Table 3 footnote: 'T1/2 of M4 in",
        "non-Asian male and female patients was 25.7 and 14.7 hours,",
        "respectively, FOR TWICE-DAILY administration', and log(2) * 1330 L /",
        "35.8 L/h = 25.7 h and log(2) * (1330 * 0.57) / 35.8 = 14.7 h using",
        "the unmodified typical values. The effect was found after the formal",
        "stepwise covariate analysis: 'After including Asian and sex,",
        "significant trends between dosing frequency (twice daily versus once",
        "daily) and EBEs of M4 were observed. The effects of dosing frequency",
        "on CLM4/FM4 and VM4/FM4 were formally tested in NONMEM and were both",
        "significant, with a P < .01'. The Discussion attributes it to",
        "'saturation in the formation/elimination of M4 under the twice-daily",
        "regimen'. NOTE the covariate therefore encodes a regimen-level",
        "empirical correction and NOT a mechanistic property; do not use this",
        "model to extrapolate to dosing intervals other than q12h and q24h."
      ),
      source_name        = "FREQ"
    ),
    OCC = list(
      description        = "Occasion index for inter-occasion variability: 1 = the first week of treatment, 2 = the remaining weeks",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "Gosselin 2015 Results 'Population PK Modeling of M4': the retained M4",
        "model carried 'IOV on all PK parameters with 2 occasions (ie, one for",
        "the first week and another for the remaining weeks)'. The motesanib",
        "model likewise retained 'IOV on CL/F' (Results 'Population PK",
        "Modeling of Motesanib', Table 3 'IOV on CL/F (%) = 18.1'), but the",
        "paper does NOT state the number or definition of the motesanib",
        "occasions; the same two-occasion week-1-versus-later split is used",
        "here, supported by the Table 4 footnote, which derives the reported",
        "exposure statistics for BOTH analytes from 'the post hoc values of",
        "first week derived with the final population PK model'. See the",
        "vignette's Assumptions and deviations section. Decomposed inside",
        "model() into the binary indicators oc1 / oc2, following the",
        "occasion-indicator expansion used by Jonsson_2011_ethambutol.R and",
        "Chen_2023_nemonoxacin.R, because rxode2 cannot simulate the",
        "'eta ~ var | occ' multi-level IOV syntax."
      ),
      source_name        = "OCC"
    )
  )

  # Covariates that Gosselin 2015 screened in the formal stepwise covariate
  # analysis but did NOT retain in either final model. They are recorded here
  # for provenance only; none is referenced in model(), and the paper reports
  # no point estimate for any of them.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at baseline",
      units       = "years",
      type        = "continuous",
      notes       = "Tested on motesanib CL/F and on M4 CLM4/FM4 (a 'small negative trend' was seen against both sets of EBEs) but not retained at P < .01. Median 59.0 y (5th-95th percentiles 38.0-75.0), Gosselin 2015 Table 2."
    ),
    CREAT = list(
      description = "Serum creatinine at baseline",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "Tested on motesanib CL/F and Vc/F and on both M4 parameters but not retained. Median 0.8 mg/dL (5th-95th percentiles 0.5-1.3), Gosselin 2015 Table 2; 6 patients had aberrant values that were imputed to the median."
    ),
    CRCL = list(
      description = "Creatinine clearance at baseline",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Reported in the demographics table as a renal-function descriptor and screened graphically, but not retained in either final model. Median 84.0 mL/min (5th-95th percentiles 47.0-153.1), Gosselin 2015 Table 2."
    ),
    BILI = list(
      description = "Total bilirubin at baseline",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Tested on motesanib CL/F and on M4 CLM4/FM4 (a 'small negative trend' against motesanib CL/F EBEs) but not retained. Median 8.55 umol/L (5th-95th percentiles 4.0-18.8), Gosselin 2015 Table 2; 6 patients had aberrant values that were imputed to the median."
    ),
    AST = list(
      description = "Aspartate aminotransferase at baseline",
      units       = "U/L",
      type        = "continuous",
      notes       = "Tested on motesanib CL/F and on M4 CLM4/FM4 but not retained. Median 22.0 U/L (5th-95th percentiles 12.5-55.0), Gosselin 2015 Table 2."
    ),
    ALT = list(
      description = "Alanine aminotransferase at baseline",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened graphically against motesanib CL/F EBEs ('small negative trends were observed between EBEs of CL/F and age, alkaline phosphatase (ALP), alanine aminotransferase, aspartate aminotransferase (AST), and total bilirubin') but not carried into the formal NONMEM covariate analysis. Median 20.0 U/L (5th-95th percentiles 9.0-56.5), Gosselin 2015 Table 2."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 451L,
    n_studies      = 8L,
    age_range      = "5th-95th percentiles 38.0-75.0 years (median 59.0)",
    weight_range   = "5th-95th percentiles 45.0-103.8 kg (median 67.41)",
    sex_female_pct = 48.6,
    race_ethnicity = c(White = 61.9, Black = 3.3, Hispanic = 2.4, Asian = 31.5, Other = 0.9),
    disease_state  = paste(
      "Patients with advanced solid tumors: 246 (54.5%) non-small cell lung",
      "cancer, 49 (10.9%) breast cancer, 43 (9.5%) gastrointestinal stromal",
      "tumor, and the remainder thyroid, pancreatic, colorectal, sarcoma,",
      "kidney, neuroendocrine and other tumor types."
    ),
    dose_range     = paste(
      "Oral motesanib 25 to 175 mg, once daily (50, 75, 100, 125 or 175 mg",
      "q24h) or twice daily (25 or 75 mg q12h), given in 28-day cycles with",
      "treatment holidays of 2 days or 1 week; plus a single oral 125 mg",
      "[14C]-labelled dose in a mass-balance study."
    ),
    regions        = "Not stated in the publication; 102 of the 142 Asian patients (72%) were Japanese",
    notes          = paste(
      "Baseline demographics from Gosselin 2015 Table 2 (n = 451). The 451",
      "patients are the union of the two analysis datasets: 445 patients",
      "contributed motesanib concentrations and 249 contributed M4",
      "concentrations, of whom 6 had no detectable motesanib concentration",
      "and so appear only in the M4 analysis. M4 was measured in 4 of the 8",
      "studies (20050200, 20080639, 20060136, 20050201). 4864 non-BLQ",
      "motesanib concentrations (2782 intensive, 2082 sparse; 95.5% of those",
      "collected) and 1108 non-BLQ M4 concentrations (300 intensive, 808",
      "sparse; 86.6%) entered the fits. BLQ values were treated as missing",
      "(202 motesanib and 36 M4 records, < 5% of samples). Estimation was",
      "FOCE with INTERACTION in NONMEM VII, with the M4 model fitted",
      "SEQUENTIALLY using the post-hoc motesanib parameters as the input",
      "driving M4 formation. Model stability was assessed with a 1000-sample",
      "bootstrap (550 and 829 successful minimisations for motesanib and M4",
      "respectively); bootstrap medians were within 3.1% of the typical",
      "values and are shown in parentheses in Table 3."
    )
  )

  ini({
    # -----------------------------------------------------------------------
    # All values are the FINAL population estimates from Gosselin 2015
    # Table 3, column "Typical Values (Medians of Bootstrap)". Bootstrap
    # medians are quoted in the trailing comments for cross-reference but are
    # NOT the encoded values.
    #
    # Two Table 3 cells are affected by the PDF's symbol font: the leading
    # multiplication sign on every categorical multiplicative factor, and the
    # leading minus sign on the ALP exponent, both render as whitespace in a
    # text extraction. The multiplicative reading is fixed by the paper's own
    # supplement ("The estimate for continuous covariate is the exponent of
    # the power function and the estimate for categorical covariate is the
    # multiplicative factor", Supplement Table 1 note, which also prints the
    # surviving glyphs "x 0.327" and "x1.52"), and the ALP sign is fixed by
    # the bootstrap median printed next to it, "(-0.218)". Both readings are
    # confirmed arithmetically against the published half-lives; see the
    # vignette's source-trace section.
    # -----------------------------------------------------------------------

    # -- Motesanib absorption: first order with a lag time, both modified by a
    #    high-fat meal. The absorption rate constant is very fast in the
    #    fasted state, consistent with the Table 4 observation that "maximum
    #    concentrations of motesanib are rapidly reached (about 30 minutes
    #    after the dose)".
    ltlag <- log(0.229)  ; label("Motesanib absorption lag time in the fasted state (h, log scale)")             # Gosselin 2015 Table 3 'Lag time (h) = 0.229' (RSE 0.3%; bootstrap median 0.229)
    lka   <- log(9.84)   ; label("Motesanib first-order absorption rate constant in the fasted state (1/h, log scale)") # Gosselin 2015 Table 3 'Ka (h-1) = 9.84' (RSE 18.3%; bootstrap median 9.97)

    # -- Motesanib disposition. All four are APPARENT values (divided by the
    #    oral bioavailability F, which the paper never separates out). CL/F is
    #    the typical value for a MAN at the reference albumin of 39 g/L; Vc/F
    #    is the typical value at ALB = 39 g/L, ALP = 101 U/L and WT = 68.3 kg.
    lcl   <- log(62.7)   ; label("Motesanib apparent clearance CL/F in men at ALB = 39 g/L (L/h, log scale)")    # Gosselin 2015 Table 3 'CL/F (L/h) = 62.7 for men' (RSE 3.5%; bootstrap median 62.7)
    lvc   <- log(240)    ; label("Motesanib apparent central volume Vc/F at ALB = 39 g/L, ALP = 101 U/L, WT = 68.3 kg (L, log scale)") # Gosselin 2015 Table 3 'Vc/F (L) = 240' (RSE 6.4%; bootstrap median 242)
    lq    <- log(93.5)   ; label("Motesanib apparent intercompartmental clearance Q/F (L/h, log scale)")         # Gosselin 2015 Table 3 'Q/F (L/h) = 93.5' (RSE 16.3%; bootstrap median 90.7)
    lvp   <- log(195)    ; label("Motesanib apparent peripheral volume Vp/F (L, log scale)")                     # Gosselin 2015 Table 3 'Vp/F (L) = 195' (RSE 6.1%; bootstrap median 192)

    # -- M4 disposition, one compartment. Both are APPARENT values divided by
    #    FM4, which the paper defines as representing "the fraction of oral
    #    absorption of motesanib, the fraction of biotransformation from
    #    motesanib to M4, and information pertaining to the molecular weight
    #    of the parent and metabolite" (Methods, 'PK Model Buildup of M4').
    #    CLM4/FM4 is the typical value in a NON-ASIAN patient on TWICE-daily
    #    dosing; VM4/FM4 is the typical value in a MAN on TWICE-daily dosing.
    lcl_m4 <- log(35.8)  ; label("M4 apparent clearance CLM4/FM4 in non-Asian patients on twice-daily dosing (L/h, log scale)") # Gosselin 2015 Table 3 'CLM4/FM4 (L/h) = 35.8 for non-Asian' (RSE 7.2%; bootstrap median 35.8)
    lvc_m4 <- log(1330)  ; label("M4 apparent volume of distribution VM4/FM4 in men on twice-daily dosing (L, log scale)")      # Gosselin 2015 Table 3 'VM4/FM4 (L) = 1330 for men' (RSE 15.1%; bootstrap median 1331)

    # -- Covariate effects on motesanib. Continuous covariates enter as power
    #    functions centred at the analysis-set medians; categorical covariates
    #    enter as multiplicative factors (Methods, 'Population PK Modeling of
    #    Motesanib').
    e_alb_cl          <- 0.971  ; label("Power exponent on (ALB/39) for motesanib CL/F (unitless)")              # Gosselin 2015 Table 3 'ALB on CL/F = (ALB/39)^0.971' (RSE 21.8%; bootstrap median 0.979)
    e_sexf_cl         <- 0.882  ; label("Multiplicative factor on motesanib CL/F for women vs men (unitless)")   # Gosselin 2015 Table 3 'Sex on CL/F = x 0.882 for women' (RSE 4.3%); 62.7 * 0.882 = 55.3 L/h, matching the Discussion's '62.7 vs 55.3 L/h'
    e_alb_vc          <- 1.66   ; label("Power exponent on (ALB/39) for motesanib Vc/F (unitless)")              # Gosselin 2015 Table 3 'ALB on Vc/F = (ALB/39)^1.66' (RSE 25.9%; bootstrap median 1.62)
    e_alp_vc          <- -0.217 ; label("Power exponent on (ALP/101) for motesanib Vc/F (unitless)")             # Gosselin 2015 Table 3 'ALP on Vc/F = (ALP/101)^-0.217' (RSE 35.9%; bootstrap median -0.218)
    e_wt_vc           <- 0.612  ; label("Power exponent on (WT/68.3) for motesanib Vc/F (unitless)")             # Gosselin 2015 Table 3 'WT on Vc/F = (WT/68.3)^0.612' (RSE 24.8%; bootstrap median 0.612)
    e_fed_highfat_ka  <- 0.0245 ; label("Multiplicative factor on motesanib Ka for a high-fat meal vs fasted (unitless)")       # Gosselin 2015 Table 3 'Diet on Ka = x 0.0245' (RSE 15.2%; bootstrap median 0.0241); 9.84 * 0.0245 = 0.241 1/h when fed
    e_fed_highfat_tlag <- 1.079 ; label("Multiplicative factor on the motesanib lag time for a high-fat meal vs fasted (unitless)") # Gosselin 2015 Table 3 'Diet on lag time = x 1.079' (RSE 7.2%; bootstrap median 1.081); 0.229 * 1.079 = 0.247 h when fed

    # -- Covariate effects on M4. All four are multiplicative factors on the
    #    reference (non-Asian, male, twice-daily) typical values.
    e_race_asian_cl_m4 <- 1.4   ; label("Multiplicative factor on M4 CLM4/FM4 for Asian vs non-Asian patients (unitless)")      # Gosselin 2015 Table 3 'Asian on CLM4/FM4 = x 1.4 for Asian' (RSE 8.4%; bootstrap median 1.40)
    e_regi_qd_cl_m4    <- 1.52  ; label("Multiplicative factor on M4 CLM4/FM4 for once-daily vs twice-daily dosing (unitless)") # Gosselin 2015 Table 3 'Dosing interval on CLM4/FM4 = x 1.52 for once-daily dose' (RSE 8.0%; bootstrap median 1.51)
    e_sexf_vc_m4       <- 0.57  ; label("Multiplicative factor on M4 VM4/FM4 for women vs men (unitless)")                      # Gosselin 2015 Table 3 'Sex on VM4/FM4 = x 0.57 for women' (RSE 12.2%; bootstrap median 0.57)
    e_regi_qd_vc_m4    <- 0.45  ; label("Multiplicative factor on M4 VM4/FM4 for once-daily vs twice-daily dosing (unitless)")  # Gosselin 2015 Table 3 'Dosing interval on VM4/FM4 = x 0.45 for once-daily dose' (RSE 15.7%; bootstrap median 0.45)

    # -----------------------------------------------------------------------
    # Inter-individual variability. Gosselin 2015 Methods: "The
    # interindividual variability (IIV) and interoccasion variability (IOV)
    # were modeled as an exponential function to positively constrain the
    # individual parameter values", so each eta is log-normal and Table 3
    # reports its magnitude as a percentage.
    #
    # SCALE OF THE REPORTED PERCENTAGES. Table 3's "IIV on <param> (%)" rows
    # are read here as the standard-deviation-scale CV, 100 * sqrt(omega^2),
    # so every variance below is (IIV% / 100)^2. Keep separate what this
    # reading PROVES from what it ASSUMES.
    #
    #   PROVEN -- the percentages are NOT on the variance scale. For a
    #      variance estimated from N subjects the asymptotic relative
    #      standard error cannot fall below sqrt(2/N); with N = 445 that
    #      floor is 6.70%. Table 3 reports 5.66% for "IIV on CL/F", below
    #      the floor for a variance but comfortably above the SD-scale floor
    #      sqrt(1/(2N)) = 3.35%. Reading 0.359 as omega^2 is therefore
    #      impossible, and omega^2 = (IIV%/100)^2 is the right shape.
    #
    #   ASSUMED -- that the plain 100 * sqrt(omega^2) is meant rather than
    #      the exact log-normal 100 * sqrt(exp(omega^2) - 1). Both are
    #      "SD-ish" and the RSE argument above does NOT separate them; the
    #      paper never prints its formula. The plain form is taken because it
    #      is the usual NONMEM / PsN reporting convention, and because it is
    #      the scale on which the authors do their own arithmetic: the
    #      Discussion's claim that the Vc/F covariates "explained 14% of
    #      total IIV on Vc/F" is exactly (57.0 - 49.0) / 57.0 from the
    #      supplement's step table, i.e. printed percentages subtracted
    #      directly.
    #
    #   IMMATERIAL EXCEPT FOR Ka. Below about 50% the two readings differ by
    #      under 4% in omega (CL/F 0.359 vs 0.346; Vc/F 0.495 vs 0.470). The
    #      choice bites only on "IIV on Ka (%) = 178", where the plain
    #      reading gives omega = 1.78 and the exact one omega = 1.195. No
    #      published statistic isolates the Ka variance -- Table 4's Cmax CV%
    #      is a post-hoc EBE spread and Ka shrinkage is not reported -- so
    #      that single variance carries genuine residual uncertainty. See the
    #      vignette's Assumptions and deviations section.
    #
    #   CORROBORATION where the choice is nearly moot: AUC(0-24) at steady
    #      state equals dose / (CL/F), so its dispersion is driven by
    #      omega_IIV,CL and omega_IOV,CL. On the encoded scale the predicted
    #      coefficient of variation is sqrt(exp(0.359^2 + 0.181^2) - 1) =
    #      41.9%, against the published 42.3% for 125 mg once daily. (Table 4
    #      values are post-hoc EBE spreads, which shrinkage deflates
    #      slightly, so a predicted CV at or just below the published one is
    #      the expected relationship.)
    #
    # NO IIV is reported on the lag time or on any covariate effect.
    # -----------------------------------------------------------------------
    etalka ~ 3.1684    # var = 1.78^2  ; Gosselin 2015 Table 3 'IIV on Ka (%) = 178' (RSE 10.4%; bootstrap median 176)
    etalcl ~ 0.128881  # var = 0.359^2 ; Gosselin 2015 Table 3 'IIV on CL/F (%) = 35.9' (RSE 5.66%; bootstrap median 35.9)
    etalvc ~ 0.245025  # var = 0.495^2 ; Gosselin 2015 Table 3 'IIV on Vc/F (%) = 49.5' (RSE 7.71%; bootstrap median 48.6)
    etalq  ~ 0.142884  # var = 0.378^2 ; Gosselin 2015 Table 3 'IIV on Q/F (%) = 37.8' (RSE 21.5%; bootstrap median 36.7)
    etalvp ~ 0.110889  # var = 0.333^2 ; Gosselin 2015 Table 3 'IIV on Vp/F (%) = 33.3' (RSE 9.82%; bootstrap median 33.1)

    # M4 IIV. Gosselin 2015 Results also reports "a correlation between
    # CLM4/FM4 and VM4/FM4" in the retained M4 model, but Table 3 does not
    # print the correlation coefficient and no supplement reports it, so the
    # two etas are left UNCORRELATED here (equivalent to fixing the
    # covariance at zero). See the vignette's Assumptions and deviations
    # section; no variance was invented.
    etalcl_m4 ~ 0.123904  # var = 0.352^2 ; Gosselin 2015 Table 3 'IIV on CLM4/FM4 (%) = 35.2' (RSE 14.1%; bootstrap median 34.7)
    etalvc_m4 ~ 0.136900  # var = 0.370^2 ; Gosselin 2015 Table 3 'IIV on VM4/FM4 (%) = 37.0' (RSE 36.0%; bootstrap median 35.7)

    # -----------------------------------------------------------------------
    # Inter-occasion variability, two occasions (week 1 versus the remaining
    # weeks). Written as the occasion-indicator expansion rather than the
    # `eta ~ var | OCC` multi-level syntax because rxode2 parses that syntax
    # but cannot SIMULATE from it (it emits an unbound THETA), and every
    # model in the library must solve. The expansion is also the more literal
    # translation of a NONMEM $PK IOV block: the first occasion's variance is
    # the estimated one and each further occasion repeats it, exactly as
    # $OMEGA BLOCK(1) SAME does. Founding precedents in this library:
    # Jonsson_2011_ethambutol.R, Aregbe_2012_alvespimycin.R,
    # Chen_2023_nemonoxacin.R, Xie_2019_agomelatine.R.
    # -----------------------------------------------------------------------
    etaiov_cl_1    ~ 0.032761         # var = 0.181^2 ; Gosselin 2015 Table 3 'IOV on CL/F (%) = 18.1' (RSE 6.72%; bootstrap median 18.1)
    etaiov_cl_2    ~ fixed(0.032761)  # same variance on occasion 2 (NONMEM $OMEGA BLOCK(1) SAME)
    etaiov_cl_m4_1 ~ 0.060025         # var = 0.245^2 ; Gosselin 2015 Table 3 'IOV on CLM4/FM4 (%) = 24.5' (RSE 20.4%; bootstrap median 24.2)
    etaiov_cl_m4_2 ~ fixed(0.060025)  # same variance on occasion 2
    etaiov_vc_m4_1 ~ 0.290521         # var = 0.539^2 ; Gosselin 2015 Table 3 'IOV on VM4/FM4 (%) = 53.9' (RSE 19.5%; bootstrap median 52.2)
    etaiov_vc_m4_2 ~ fixed(0.290521)  # same variance on occasion 2

    # -----------------------------------------------------------------------
    # Residual variability. Gosselin 2015 Methods: "Residual variability was
    # described using a statistical model with an additive and a proportional
    # component." Both components of both error models were held FIXED during
    # estimation -- Table 3 prints the literal token "FIX" in the relative
    # standard error column for all four rows, in place of a percentage.
    # The additive terms are in the paper's concentration unit, ug/L.
    # -----------------------------------------------------------------------
    addSd     <- fixed(0.334) ; label("Motesanib additive residual error (SD, ug/L)")        # Gosselin 2015 Table 3 motesanib 'Additive error (ug/L) = 0.334', RSE column 'FIX'
    propSd    <- fixed(0.473) ; label("Motesanib proportional residual error (SD, fraction)") # Gosselin 2015 Table 3 motesanib 'Proportional error (%) = 47.3', RSE column 'FIX'
    addSd_m4  <- fixed(4.82)  ; label("M4 additive residual error (SD, ug/L)")               # Gosselin 2015 Table 3 M4 'Additive error (ug/L) = 4.82', RSE column 'FIX'
    propSd_m4 <- fixed(0.448) ; label("M4 proportional residual error (SD, fraction)")        # Gosselin 2015 Table 3 M4 'Proportional error (%) = 44.8%', RSE column 'FIX'
  })

  model({
    # 1. Unit bridge. Doses are in mg and volumes in L, so an amount over a
    #    volume is mg/L; the paper's concentrations and additive residual
    #    errors are in ug/L.
    ugPerMg <- 1000

    # 2. Occasion indicators for the two-occasion IOV (1 = first week of
    #    treatment, 2 = remaining weeks).
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    iov_cl    <- oc1 * etaiov_cl_1    + oc2 * etaiov_cl_2
    iov_cl_m4 <- oc1 * etaiov_cl_m4_1 + oc2 * etaiov_cl_m4_2
    iov_vc_m4 <- oc1 * etaiov_vc_m4_1 + oc2 * etaiov_vc_m4_2

    # 3. Individual motesanib absorption parameters. The high-fat meal
    #    multiplies both, so FED_HIGHFAT = 0 recovers the fasted typical
    #    values exactly.
    ka   <- exp(lka + etalka) * e_fed_highfat_ka^FED_HIGHFAT
    tlag <- exp(ltlag)        * e_fed_highfat_tlag^FED_HIGHFAT

    # 4. Individual motesanib disposition parameters. Continuous covariates
    #    are power functions centred at the analysis-set medians; the sex
    #    effect is a multiplicative factor relative to men.
    cl <- exp(lcl + etalcl + iov_cl) * (ALB / 39)^e_alb_cl * e_sexf_cl^SEXF
    vc <- exp(lvc + etalvc) * (ALB / 39)^e_alb_vc * (ALP / 101)^e_alp_vc * (WT / 68.3)^e_wt_vc
    q  <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)

    # 5. Individual M4 disposition parameters, relative to the reference
    #    non-Asian man on twice-daily dosing.
    cl_m4 <- exp(lcl_m4 + etalcl_m4 + iov_cl_m4) *
      e_race_asian_cl_m4^RACE_ASIAN * e_regi_qd_cl_m4^REGI_QD
    vc_m4 <- exp(lvc_m4 + etalvc_m4 + iov_vc_m4) *
      e_sexf_vc_m4^SEXF * e_regi_qd_vc_m4^REGI_QD

    # 6. ODE system. The M4 model was fitted SEQUENTIALLY, with the post-hoc
    #    motesanib profile as the input driving M4 formation (Methods, 'PK
    #    Model Buildup of M4'), so the formation flux into central_m4 is the
    #    WHOLE apparent elimination flux out of motesanib's central
    #    compartment. That is what makes CLM4 and VM4 apparent values scaled
    #    by FM4 = F * fm * (MW ratio): with the parent state holding
    #    amount/F and the metabolite state holding amount/FM4, the transfer
    #    term is exactly cl * central / vc with no explicit fm or
    #    molecular-weight factor. The steady-state consequence, AUC_M4 =
    #    dose / CLM4/FM4 independent of the parent clearance, reproduces
    #    Table 4: 125 mg / 54.4 L/h (once-daily non-Asian) = 2.30 mg*h/L
    #    against a published mean of 2.27 mg*h/L.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <- ka * depot -
                           cl * central / vc -
                           q * central / vc + q * peripheral1 / vp
    d/dt(peripheral1) <- q * central / vc - q * peripheral1 / vp
    d/dt(central_m4)  <- cl * central / vc - cl_m4 * central_m4 / vc_m4

    # 7. Absorption lag on the depot.
    alag(depot) <- tlag

    # 8. Observations, both in ug/L.
    Cc    <- ugPerMg * central    / vc
    Cc_m4 <- ugPerMg * central_m4 / vc_m4

    Cc    ~ add(addSd)    + prop(propSd)
    Cc_m4 ~ add(addSd_m4) + prop(propSd_m4)
  })
}
