Hopkins_2015_leflunomide <- function() {
  description <- "Parametric time-to-event (TTE) model for cessation of daily oral leflunomide due to toxicity in adult rheumatoid-arthritis patients. Two piecewise-exponential hazards describe (i) the toxicity-cessation hazard h0(t) with four intervals (< 50, 50-112, 112-204, > 204 days after leflunomide initiation) and (ii) a parallel random-censoring hazard h0ran(t) with five intervals (< 147, 147-210, 210-314, 314-350, > 350 days). CYP1A2 rs762551 C-allele carrier status is the single retained covariate: it multiplies the cessation hazard by (1 + 1.29 * SNP_CYP1A2_RS762551_C_CARRIER), a 2.29-fold increase for C-carriers vs AA homozygotes. Teriflunomide steady-state trough concentrations (total and free), leflunomide dose, CYP2C19 phenotype, and 20 other covariates were screened but not retained in the final model. Outputs sur / hazard / cumhaz (cessation) and sur_cens / hazard_cens / cumhaz_cens (random censoring) are exposed for forward simulation."
  reference <- paste(
    "Hopkins AM, Wiese MD, Proudman SM, O'Doherty CE, Upton RN, Foster DJR.",
    "Genetic polymorphism of CYP1A2 but not total or free teriflunomide",
    "concentrations is associated with leflunomide cessation in rheumatoid",
    "arthritis.",
    "Br J Clin Pharmacol. 2016;81(1):113-123 (Accepted Article published online",
    "29 August 2015).",
    "doi:10.1111/bcp.12760.",
    sep = " "
  )
  vignette <- "Hopkins_2015_leflunomide"
  units <- list(
    time          = "day",
    dosing        = "n/a (no drug-dosing events; the time-to-event model has no PK/PD parameters and does not carry the leflunomide / teriflunomide PK from the upstream Hopkins semi-physiological model)",
    concentration = "probability (the model outputs sur and sur_cens are survival probabilities, not drug concentrations)"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    cumhaz      = list(analyte = "toxicity-cessation hazard h0(t)", units = NA_character_, specimen = "not applicable", verified = FALSE),
    cumhaz_cens = list(analyte = "parallel random-censoring hazard h0ran(t)", units = NA_character_, specimen = "not applicable", verified = FALSE)
  )

  covariateData <- list(
    SNP_CYP1A2_RS762551_C_CARRIER = list(
      description        = "Binary indicator: 1 = carrier of at least one C allele at CYP1A2 rs762551 (the paper's 'C163A' promoter-region SNP; genotype AC or CC), 0 = homozygous AA (the CYP1A2*1F reference / ultra-inducible form).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (homozygous AA at CYP1A2 rs762551; the CYP1A2*1F/*1F reference)",
      notes              = "Time-fixed per subject (germline genotype). AC heterozygous and CC homozygous carriers were pooled by the authors because the CC homozygous frequency (9/105 = 8.6%) was too low for a separate phenotype layer (Discussion). Distribution in the Hopkins 2015 cohort (Table 1): AA 54/105 (51.4%), AC or CC 51/105 (48.6%). Effect direction: C-carriers have a 2.29-fold increase in the instantaneous leflunomide-cessation hazard versus AA reference (Table 3, theta10 = 1.29; lambda_CYP1A2 = 1 + theta10).",
      source_name        = "CYP1A2 C163A"
    )
  )

  # All 22 covariates that the authors screened but did not retain in the
  # final model. Table 4 lists the 17 covariates with reported delta-OFVs;
  # the Table 4 footnote lists the additional 12 covariates that produced
  # no OFV improvement. Kept here for provenance only.
  # checkModelConventions() treats covariatesDataExcluded as documentation:
  # entries must not be referenced in model().
  covariatesDataExcluded <- list(
    ALT = list(
      description        = "Alanine transaminase (measured at cessation or censor date).",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened as an exponential function normalised by the median (paper Equation 8); univariate dOFV = 3.274, P = 0.07 (Table 4). Not retained.",
      source_name        = "ALT"
    ),
    AST = list(
      description        = "Aspartate transaminase (measured at cessation or censor date).",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened as exponential-normalized; univariate dOFV = 2.849, P = 0.091 (Table 4). Not retained.",
      source_name        = "AST"
    ),
    LOADING_DOSE = list(
      description        = "Use of a leflunomide loading dose (3 daily doses of 100 mg followed by 20 mg/day). 1 = loading dose used, 0 = no loading dose.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no loading dose)",
      notes              = "Univariate dOFV = 2.784, P = 0.095 (Table 4). Not retained. Only 5/105 (4.8%) subjects received a loading dose.",
      source_name        = "Use of leflunomide loading dose"
    ),
    SEXF = list(
      description        = "Patient sex indicator (1 = female, 0 = male).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Univariate dOFV = 2.242, P = 0.134 (Table 4). Not retained. 80/105 (76.2%) of the cohort was female. Encoded against the canonical SEXF = 1 (female) convention.",
      source_name        = "Gender"
    ),
    LEF_DAILY_DOSE = list(
      description        = "Average daily leflunomide dose across the study period.",
      units              = "mg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened as exponential-normalized on the median; univariate dOFV = 1.953, P = 0.162 (Table 4). Median dose 14.40 mg/day (range 6.64-100.00). Not retained.",
      source_name        = "Average daily leflunomide dose"
    ),
    DHODH_HAP2 = list(
      description        = "DHODH haplotype II carrier indicator (1 = carrier, 0 = non-carrier).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-carrier)",
      notes              = "Univariate dOFV = 1.731, P = 0.188 (Table 4). Not retained. 56/105 (53.3%) were haplotype-II carriers.",
      source_name        = "DHODH haplotype II status (dichotomous)"
    ),
    SSZ_DOSE = list(
      description        = "Concomitant daily sulfasalazine dose (measured at cessation or censor date).",
      units              = "mg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened as a linear function of dose (paper Equation 9); univariate dOFV = 1.295, P = 0.255 (Table 4). Median 2000 mg/day (range 0-3000). Not retained.",
      source_name        = "SSZ dose"
    ),
    ANTI_CCP_POS = list(
      description        = "Anti-cyclic citrullinated peptide antibody positivity at RA diagnosis (1 = positive, 0 = negative).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (negative)",
      notes              = "Univariate dOFV = 1.164, P = 0.281 (Table 4). Not retained. 60/103 (58.2%) tested positive.",
      source_name        = "Anti-CCP positivity"
    ),
    DHODH_RS3213422 = list(
      description        = "DHODH rs3213422 genotype (three-level polychotomous: AA / AC / CC).",
      units              = "(categorical)",
      type               = "categorical",
      reference_category = "AC (the most frequent category: 53/105)",
      notes              = "Polychotomous form screened; univariate dOFV = 0.862, P = 0.353, 2 df (Table 4). Not retained. Distribution: AA 24/105 (22.9%), AC 53/105 (50.5%), CC 28/105 (26.7%).",
      source_name        = "DHODH rs3213422 genotype (polychotomous)"
    ),
    CYP2C19_PHENOTYPE = list(
      description        = "CYP2C19 predicted metabolic phenotype (four-level: ultra-rapid / extensive / poor-or-intermediate / unknown). Unknown metabolizers were imputed with the mode. Determined from rs4244285 and rs12248560 genotypes per Wiese et al. 2012.",
      units              = "(categorical)",
      type               = "categorical",
      reference_category = "Extensive metabolizer (the most frequent category: 51/105)",
      notes              = "Polychotomous form screened; univariate dOFV = 1.653, P = 0.438, 3 df (Table 4). Not retained. Distribution: ultra-rapid 26/105 (24.8%), extensive 51/105 (48.6%), poor/intermediate 20/105 (19.0%), unknown 8/105 (7.6%). This finding does NOT confirm the earlier report of Wiese et al. 2012 (which used a Cox proportional-hazards model in a smaller subset of the same cohort).",
      source_name        = "CYP2C19 'phenotype' (polychotomous)"
    ),
    FREE_TERI_SS_TROUGH = list(
      description        = "Predicted free teriflunomide steady-state trough concentration after 365 daily doses at the individual's average leflunomide dose. Derived from a semi-physiological pharmacokinetic model (Hopkins 2015 reference [37]).",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened as exponential-normalized; univariate dOFV = 0.487, P = 0.485 (Table 4). Median 0.0436 mg/L (range 0.0022-0.2879). Not retained. This was the paper's headline exposure metric because free (unbound) teriflunomide is the pharmacologically active fraction (> 99% protein-bound total teriflunomide). Population parameter estimates (CLint, Vbody, fu) imputed via median (adjusted for FFM and ALT) in 36/105 subjects without blood samples.",
      source_name        = "Predicted free teriflunomide steady-state trough concentrations"
    ),
    AGE = list(
      description        = "Patient age at leflunomide initiation.",
      units              = "year",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened as exponential-normalized on the median (58.4 y); univariate dOFV = 0.473, P = 0.492 (Table 4). Not retained. Range 19.2-85.9 y.",
      source_name        = "Age"
    ),
    TRIPLE_THERAPY = list(
      description        = "1 = subject was on triple therapy (methotrexate + sulfasalazine + hydroxychloroquine) at leflunomide initiation, 0 = not on triple therapy.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not on triple therapy at initiation)",
      notes              = "Univariate dOFV = 0.388, P = 0.533 (Table 4). Not retained. 59/105 (56.2%) were on triple therapy at initiation.",
      source_name        = "Triple therapy at initiation"
    ),
    MTX_DOSE = list(
      description        = "Concomitant weekly methotrexate dose (measured at cessation or censor date).",
      units              = "mg/week",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened as a linear function of dose; univariate dOFV = 0.365, P = 0.546 (Table 4). Median 20 mg/wk (range 0-25). Not retained.",
      source_name        = "MTX dose"
    ),
    HCQ_DOSE = list(
      description        = "Concomitant daily hydroxychloroquine dose (measured at cessation or censor date).",
      units              = "mg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened as a linear function of dose; univariate dOFV = 0.326, P = 0.568 (Table 4). Median 400 mg/day (range 0-800). Not retained.",
      source_name        = "HCQ dose"
    ),
    RF_POS = list(
      description        = "Rheumatoid factor positivity at RA diagnosis (1 = positive, 0 = negative).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (negative)",
      notes              = "Univariate dOFV = 0.121, P = 0.728 (Table 4). Not retained. 67/104 (64.4%) tested positive.",
      source_name        = "RF positivity"
    ),
    NULL_COVARIATES_NO_OFV = list(
      description        = "Composite entry for the 12 covariates whose univariate screen produced no OFV improvement (paper Table 4 footnote): baseline DAS28, shared-epitope status, weight, height, fat-free mass (FFM), smoking status, ABCG2 rs2231142 genotype, creatinine clearance (CLcr, Cockcroft-Gault/IBW), fraction unbound (fu), intrinsic clearance (CLint), volume of distribution (Vbody), and predicted total teriflunomide steady-state trough concentrations. Kept here for provenance only; no per-covariate values are reported in the paper because the univariate step returned no improvement.",
      units              = "(mixed)",
      type               = "categorical",
      reference_category = NULL,
      notes              = "Paper Table 4 footnote: 'The covariates of baseline DAS28, SE status, weight, height, FFM, smoking status, ABCG2 genotype, CLcr, fu, CLint, Vbody and predicted total teriflunomide steady-state trough concentrations did not result in any improvement in the OBJ.' A backward-elimination step (P < 0.005) was planned but not necessary given only one covariate was retained.",
      source_name        = "Paper Table 4 footnote"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 105L,
    n_studies       = 1L,
    age_range       = "19.2-85.9 years (median 58.4 at leflunomide initiation)",
    weight_range    = "40.4-148.5 kg (median 72.0 at cessation or censor date)",
    sex_female_pct  = 76.2,
    race_ethnicity  = "not reported (single-centre Adelaide cohort; the paper does not disaggregate ethnicity)",
    disease_state   = "Adult rheumatoid arthritis. DMARD-naive at RA diagnosis per revised ACR Criteria; enrolled in the Early Arthritis inception cohort at the Royal Adelaide Hospital between 2000 and 2013. Leflunomide was added after triple therapy (methotrexate + sulfasalazine + hydroxychloroquine) failed to control disease activity. 60/103 (58.2%) were anti-CCP positive at diagnosis, 67/104 (64.4%) RF positive, 71/101 (70.3%) SE positive.",
    dose_range      = "Daily oral leflunomide. First 3 years of the cohort: 3 daily doses of 100 mg loading + 20 mg/day maintenance (5/105 patients only). Subsequently: 10 mg/day initial, up-titrated to 20 mg/day in persistent disease. Median daily dose across the study period 14.40 mg/day (range 6.64-100.00).",
    regions         = "Australia (Royal Adelaide Hospital; single-centre)",
    observation_window = "Up to 60 weeks after leflunomide initiation; censored at week 52 or at the closest clinic visit if > 52 weeks. 33/105 (31.4%) had their last clinic visit before day 365; 26/105 (24.8%) between days 365 and 420.",
    co_medication   = "Concomitant DMARD dosing at cessation or censor date: methotrexate (MTX) in 80/105 (76.1%; median 20 mg/wk, range 0-25), sulfasalazine (SSZ) in 70/105 (66.7%; median 2000 mg/day, range 0-3000), hydroxychloroquine (HCQ) in 91/105 (86.7%; median 400 mg/day, range 0-800).",
    notes           = "Observational retrospective cohort. 115 patients commenced leflunomide; 10 excluded for incomplete records, leaving N = 105 eligible for the time-to-event analysis. 34/105 (32.4%) ceased leflunomide due to toxicity before day 365; 3/105 (2.9%) ceased for non-toxicity reasons; 9/105 (8.6%) had another DMARD added for inadequate response; the remainder were censored at their last clinic visit. Toxicities driving cessation (Table 2): gastrointestinal 20, neutropenia 7, fatigue/dizziness/headaches 5, elevated liver enzymes 4, respiratory symptoms including pneumonitis 4, other 8. Teriflunomide (leflunomide's active metabolite) steady-state trough concentrations (total and free) were derived from a per-subject semi-physiological pharmacokinetic model (Hopkins 2015 reference [37]) and used only in the covariate screen; they are not carried in the packaged TTE model. Blood samples were available from only 69/105 subjects; the remaining 36 had per-subject PK parameter estimates imputed with the population median adjusted for FFM and ALT."
  )

  ini({
    # All hazards on the natural log scale. Estimates come from Table 3,
    # covariate-model column (theta10 = 1.29 CYP1A2 effect retained;
    # baseline hazards recorded here are the covariate-model values, NOT
    # the structural-model values from the Table 3 left half). Bootstrap
    # 95% CI and %RSE for each parameter are captured in the vignette
    # source-trace table.

    # Cessation hazard h0(t): piecewise-exponential, four time intervals.
    lh1 <- log(0.00161); label("Log cessation hazard for t < 50 days after leflunomide initiation (1/day)")            # Table 3 covariate-model theta1
    lh2 <- log(0.00034); label("Log cessation hazard for 50 <= t < 112 days after leflunomide initiation (1/day)")     # Table 3 covariate-model theta2
    lh3 <- log(0.00151); label("Log cessation hazard for 112 <= t < 204 days after leflunomide initiation (1/day)")    # Table 3 covariate-model theta3
    lh4 <- log(0.00021); label("Log cessation hazard for t >= 204 days after leflunomide initiation (1/day)")          # Table 3 covariate-model theta4

    # Random-censoring hazard h0ran(t): piecewise-exponential, five time intervals.
    lhran1 <- log(0.00024); label("Log random-censoring hazard for t < 147 days (1/day)")               # Table 3 covariate-model theta5
    lhran2 <- log(0.00142); label("Log random-censoring hazard for 147 <= t < 210 days (1/day)")        # Table 3 covariate-model theta6
    lhran3 <- log(0.00068); label("Log random-censoring hazard for 210 <= t < 314 days (1/day)")        # Table 3 covariate-model theta7
    lhran4 <- log(0.00775); label("Log random-censoring hazard for 314 <= t < 350 days (1/day)")        # Table 3 covariate-model theta8
    lhran5 <- log(0.03440); label("Log random-censoring hazard for t >= 350 days (1/day)")              # Table 3 covariate-model theta9

    # CYP1A2 covariate effect on the cessation hazard.
    # Paper parameterisation (Table 3 caption): in C-allele carriers,
    #   lambda_CYP1A2 = 1 + theta10 = 1 + 1.29 = 2.29,
    # which is applied multiplicatively on the baseline hazard, i.e. the
    # cessation hazard for C-carriers is (1 + theta10) * h0(t).
    e_cyp1a2_haz <- 1.29; label("CYP1A2 C-allele carrier additive-multiplicative hazard shift (unitless; hazard multiplier is 1 + e_cyp1a2_haz * SNP_CYP1A2_RS762551_C_CARRIER)")  # Table 3 covariate-model theta10

    # No IIV. Paper Methods: "As time-to-event models only use one
    # observation for each individual, random effects on the baseline
    # hazard could not be estimated, i.e. the same baseline hazard is
    # assumed for all subjects." No eta* parameters are added.

    # No residual error. NONMEM run uses the parametric survival
    # likelihood via the F_FLAG / LIKE mechanism (paper Methods
    # Equation 3 defines the probability density f(t) = S(t) * h(t)).
    # This nlmixr2 translation is intended for forward simulation;
    # hazard, cumhaz, sur and their _cens twins are exposed as
    # derived outputs.
  })
  model({
    # Back-transformed piecewise cessation hazards (1/day).
    h1 <- exp(lh1)
    h2 <- exp(lh2)
    h3 <- exp(lh3)
    h4 <- exp(lh4)

    # Back-transformed piecewise random-censoring hazards (1/day).
    hran1 <- exp(lhran1)
    hran2 <- exp(lhran2)
    hran3 <- exp(lhran3)
    hran4 <- exp(lhran4)
    hran5 <- exp(lhran5)

    # Cessation baseline hazard step function with break-points at
    # 50, 112, 204 days (paper Table 3 covariate-model rows). Boundaries
    # follow the paper's "50 <= days < 112" convention (i.e., strict-less
    # comparison so the right-hand endpoint 50 falls in the NEXT interval).
    h0 <- ifelse(t < 50,  h1,
          ifelse(t < 112, h2,
          ifelse(t < 204, h3, h4)))

    # CYP1A2 covariate multiplier on the cessation hazard.
    lambda_cyp1a2 <- 1 + e_cyp1a2_haz * SNP_CYP1A2_RS762551_C_CARRIER

    # Effective cessation hazard and integrated survival state.
    hazard <- h0 * lambda_cyp1a2
    d/dt(cumhaz) <- hazard
    cumhaz(0) <- 0
    sur <- exp(-cumhaz)

    # Random-censoring hazard step function with break-points at
    # 147, 210, 314, 350 days (paper Table 3 covariate-model rows). No
    # covariates enter the censoring hazard.
    hazard_cens <- ifelse(t < 147, hran1,
                   ifelse(t < 210, hran2,
                   ifelse(t < 314, hran3,
                   ifelse(t < 350, hran4, hran5))))
    d/dt(cumhaz_cens) <- hazard_cens
    cumhaz_cens(0) <- 0
    sur_cens <- exp(-cumhaz_cens)
  })
}
