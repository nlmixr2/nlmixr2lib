Jin_2025_benralizumab <- function() {
  description <- "Two-compartment population PK model of benralizumab (anti-IL-5R alpha) with first-order subcutaneous absorption in Chinese, Asian and non-Asian adults and adolescents with severe eosinophilic asthma plus healthy volunteers (Jin 2025), pooling 12 phase I-III studies; body weight on CL/V2/V3, Asian race on CL, anti-drug antibody on CL, study- and dose-specific absolute subcutaneous bioavailability, and study-stratified log-scale residual error"
  reference <- paste(
    "Jin Y, Guiastrennec B, Stuke M, Yao Y, Zhang Y, Barker P, Jison M,",
    "Penland RC, Ding J, Lukka PB.",
    "Population pharmacokinetics and exposure-response analysis of benralizumab in",
    "Chinese adults, adolescents, and pediatric participants with severe eosinophilic asthma.",
    "Clin Pharmacokinet. 2025;64:1233-1245. doi:10.1007/s40262-025-01538-9.",
    "Parameter values are the final-model estimates in Resource 10 of the electronic",
    "supplementary material. Updates the legacy model of Yan L, Wang B, Chia YL, Roskos LK.",
    "Clin Pharmacokinet. 2019;58:943-58; doi:10.1007/s40262-019-00738-4.",
    sep = " "
  )
  vignette <- "Jin_2025_benralizumab"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot       = list(analyte = "benralizumab", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "benralizumab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "benralizumab", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power (multiplicative) effects normalised to 70 kg on all three disposition parameters:",
        "(WT/70)^0.849 on CL, (WT/70)^0.799 on V2 and (WT/70)^0.639 on V3 (Jin 2025 Resource 10,",
        "'Parameter-covariate relationships'; the 70 kg centering is stated in the Resource 10",
        "footnote, 'BWGT baseline body weight (centered around 70 kg)'). Baseline value, time-fixed",
        "per subject. Jin 2025 Table 2 gives 77.5 +/- 18.9 kg (40.3-204) in adults and",
        "60.6 +/- 20.6 kg (40-155) in adolescents. Source column BWGT."
      ),
      source_name        = "BWGT"
    ),
    RACE_ASIAN = list(
      description        = "Asian race indicator, 1 = Asian (all self-reported Asian participants, including all Chinese participants), 0 = non-Asian",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Asian)",
      notes              = paste(
        "Log-scale multiplicative effect on clearance, CL * exp(0.0952 * RACE_ASIAN), i.e. a 9.99%",
        "higher CL in Asian participants (Jin 2025 Resource 10 row 'Beta_CL, ASIAN_1'; the +9.99%",
        "transformed value follows the Resource 10 footnote a convention",
        "TVALUE = 100 * (exp(VALUE) - 1)). Jin 2025 Methods 2.1: 'the covariate Asian was based on the",
        "subject's race. All Chinese subjects were classified as Asian.' 590 of 2855 participants",
        "(20.7%) were Asian. Jin 2025 retained this effect even though it does not meet the",
        "clinical-relevance threshold (<10%), because it improved predictive performance for the",
        "MIRACLE study; the legacy Yan 2019 model omitted it for the same <10% reason.",
        "Jin 2025 also fitted a CHINESE-on-CL variant (beta 0.115, Resource 11 third column) but",
        "selected the ASIAN variant as final -- see covariatesDataExcluded$RACE_CHINESE.",
        "Source column ASIAN."
      ),
      source_name        = "ASIAN"
    ),
    ADA_POS = list(
      description        = "Anti-drug antibody status, 1 = positive ADA titer, 0 = no ADA",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no anti-drug antibody detected)",
      notes              = paste(
        "Log-scale multiplicative effect on clearance, CL * exp(0.762 * ADA_POS), i.e. a 114% higher",
        "CL (2.14-fold) in ADA-positive participants (Jin 2025 Resource 10 row 'RCLADA', transformed",
        "value +114% under the footnote a convention). Carries its own small NORMALLY distributed",
        "(additive, not lognormal) inter-individual random effect on the coefficient, omega 0.036 fixed",
        "-- see the etae_ada_cl comment in ini(). Time-varying in the source dataset (ADA is assessed",
        "at each predesignated visit). Jin 2025 Resource 1 defines the column as",
        "'ADA Anti-drug antibody status (0: no ADA, 1: positive ADA titer)'. The sibling",
        "modellib('Wang_2017_benralizumab') uses a HIGH-TITER (>=400) ADA definition instead;",
        "the two are not interchangeable. Source column ADA."
      ),
      source_name        = "ADA"
    ),
    STUDY_MICP220 = list(
      description        = "Study MI-CP220 indicator, 1 = participant enrolled in the phase II MI-CP220 benralizumab study, 0 = any other study",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (any of the other 11 pooled studies)",
      notes              = paste(
        "Selects a study-specific absolute subcutaneous bioavailability of 0.457 in place of the",
        "reference 0.539 (Jin 2025 Resource 10 row 'Fa1S220 (fraction)'; the Resource 10 footnote",
        "defines 'Fa1S220 absolute SC bioavailability for study 220'). Fa1S220 carries its own",
        "lognormal IIV (omega 0.373) distinct from the reference Fa1 IIV. MI-CP220 is also the single",
        "study with its own residual-error magnitude (Error_ADD2 = 0.549 log(ng/mL) versus 0.175 for",
        "the other early studies), encoded here as the separate Cc_micp220 endpoint. Time-fixed per",
        "subject. Derived from the source STUDYN column."
      ),
      source_name        = "STUDYN"
    ),
    STUDY_AMES = list(
      description        = "AMES study indicator, 1 = participant enrolled in the phase I AMES autoinjector-versus-prefilled-syringe study (NCT02968914), 0 = any other study",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (any of the other 11 pooled studies)",
      notes              = paste(
        "Selects a study-specific absolute subcutaneous bioavailability of 0.688 in place of the",
        "reference 0.539 (Jin 2025 Resource 10 row 'FaslS30 (fraction)'; the Resource 10 footnote",
        "defines 'Fa1S30 absolute SC bioavailability for AMES study'). Fa1S30 carries its own",
        "lognormal IIV (omega 0.144). This effect exists because external validation of the legacy",
        "model under-predicted the AMES concentrations: Jin 2025 Results 3.2 evaluated injection site,",
        "healthy status and an AMES study effect on Fa1, CL and V2, and 'the model featuring an effect",
        "of AMES on Fa1 performed best and was thus selected as the base model for the covariate",
        "modeling step'. AMES enrolled 180 healthy volunteers receiving a single 30 mg SC dose",
        "(Jin 2025 Table 1). Time-fixed per subject. Derived from the source STUDYN column."
      ),
      source_name        = "STUDYN"
    ),
    DOSE_HIGH = list(
      description        = "Highest-dose-cohort indicator, 1 = participant received the 200 mg subcutaneous benralizumab dose, 0 = any lower dose",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (any subcutaneous dose below 200 mg)",
      notes              = paste(
        "Step-function reduction of absolute subcutaneous bioavailability at the top of the dose range:",
        "Fa1 * exp(-0.554 * DOSE_HIGH), a 42.5% reduction (Jin 2025 Resource 10 row 'Rfa1Dose',",
        "transformed value -42.5% under the footnote a convention; the Resource 10 footnote defines",
        "'RFa1Dose relative effect of dose 200 on Fa1'). 200 mg is the highest subcutaneous dose in",
        "the pooled dataset -- Jin 2025 Methods 2.2 records that two phase II studies administered",
        "2-200 mg subcutaneously -- so DOSE_HIGH = 1 identifies the 200 mg arms only and every other",
        "cohort (2-100 mg SC, and the 0.0003-3 mg/kg intravenous phase I/II cohorts) sits in the",
        "reference group. Carries a small NORMALLY distributed (additive) random effect on the",
        "coefficient, omega 0.029 fixed -- see the etae_dosehigh_fdepot comment in ini(). Same",
        "step-at-the-top-of-the-dose-range semantics on bioavailability as",
        "modellib('Maleki_2024_brepocitinib') and modellib('Hughes_2022_brepocitinib'), except that",
        "the step here is downward. Time-fixed per subject. Derived from the assigned dose level."
      ),
      source_name        = "DOSE"
    )
  )

  # Covariates that Jin 2025 (or the legacy Yan 2019 analysis it updates) screened
  # but did not retain in the final model, plus the one covariate that was fitted
  # in a rejected model variant. Documented for provenance only; deliberately
  # absent from model().
  covariatesDataExcluded <- list(
    STUDY_PHASE3 = list(
      description = "Phase III study indicator, 1 = participant enrolled in SIROCCO, CALIMA, ZONDA, BISE or MIRACLE, 0 = participant enrolled in one of the phase I/II studies",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "DATASET-CONSTRUCTION GUIDANCE, not a model covariate. It records which of the three",
        "study-stratified residual-error magnitudes applies to a concentration record, and therefore",
        "which of the three model endpoints the record belongs to. Jin 2025 Resource 10 reports three",
        "additive-on-the-log-scale residual errors and its footnote assigns them:",
        "'ADD1 additive error (all observations except phase III and study MI-CP220),",
        "ADD2 additive error (study MI-CP220), ADD3 additive error (SIROCCO, CALIMA, ZONDA, BISE, and",
        "MIRACLE studies)'. Records with STUDY_PHASE3 = 1 are the Cc endpoint (expSd = 0.367); records",
        "with STUDY_MICP220 = 1 are the Cc_micp220 endpoint (expSd_Cc_micp220 = 0.549); all remaining",
        "records are the Cc_early endpoint (expSd_Cc_early = 0.175). It is deliberately NOT referenced",
        "in model(), because nlmixr2 selects a residual error by endpoint (cmt/dvid) rather than by a",
        "covariate value; the stratum is expressed by which endpoint a row is assigned to. Time-fixed",
        "per subject. Derived from the source STUDYN column."
      )
    ),
    RACE_CHINESE = list(
      description = "Chinese-heritage race indicator (all participants from mainland China and Taiwan)",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Fitted as an alternative to RACE_ASIAN on clearance (beta 0.115, RSE 20%, Jin 2025",
        "Resource 11 third column, giving exp(0.115) - 1 = +12.2% CL) but NOT retained in the final",
        "popPK model. Jin 2025 Results 3.2: covariate-selection results 'were similar between the",
        "models built based on Asian and Chinese participants. In both cases, an effect on CL was the",
        "most significant, with a drop in objective function value of 43 and 35 points for Asian and",
        "Chinese participants, respectively', and the Discussion adds that 'evaluation of the Chinese",
        "race covariate on top of the selected Asian race covariate effect on CL did not significantly",
        "improve the model'. The two covariates are highly correlated (Phi = 64.2%) because all",
        "Chinese participants are Asian by definition. The final model uses RACE_ASIAN.",
        "RACE_CHINESE IS retained in the exposure-response layer -- see",
        "modellib('Jin_2025_benralizumab_aaer'), where it scales the asthma-exacerbation Emax."
      )
    ),
    AGE = list(
      description = "Baseline subject age",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Screened (Jin 2025 Resource 1 and Resource 2; the legacy covariate set included age and an",
        "adult/adolescent age-group flag) and not retained. Jin 2025 Discussion reports 'the absence",
        "of age-dependent differences in pharmacokinetics'; the adults-only model extrapolated to",
        "adolescents with good predictive performance (Jin 2025 Figure 3), which is the external",
        "validation that justified omitting an age term. Body weight carries the paediatric /",
        "adolescent size effect instead."
      )
    ),
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened on CL and V2 (Jin 2025 Resource 1, column SEXF) and not retained in the final model."
    ),
    ALB = list(
      description = "Baseline serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened (Jin 2025 Resource 1 / Resource 2 legacy covariate list) and not retained."
    ),
    CRCL = list(
      description = "Baseline creatinine clearance (Cockcroft-Gault)",
      units       = "mL/min",
      type        = "continuous",
      notes       = paste(
        "Screened (Jin 2025 Resource 1, column CRCL) and not retained -- expected for a 150 kDa IgG1",
        "monoclonal antibody, which is not renally eliminated. Jin 2025 Table 2 reports",
        "109 +/- 34 mL/min (2.41-349) in adults."
      )
    ),
    ALT = list(
      description = "Baseline alanine aminotransferase",
      units       = "ukat/L",
      type        = "continuous",
      notes       = "Hepatic marker screened (Jin 2025 Resource 1 / Resource 2) and not retained."
    ),
    AST = list(
      description = "Baseline aspartate aminotransferase",
      units       = "ukat/L",
      type        = "continuous",
      notes       = "Hepatic marker screened (Jin 2025 Resource 1 / Resource 2) and not retained."
    ),
    TBIL = list(
      description = "Baseline total bilirubin",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Hepatic marker screened (Jin 2025 Resource 1, column TBL) and not retained."
    ),
    EOS = list(
      description = "Baseline blood eosinophil count",
      units       = "cells/uL",
      type        = "continuous",
      notes       = paste(
        "Screened on the PK parameters (Jin 2025 Resource 1, column BEOSL) and not retained. Baseline",
        "eosinophil count IS a retained covariate in the FEV1 exposure-response layer -- see",
        "modellib('Jin_2025_benralizumab_fev1')."
      )
    ),
    SMOKE = list(
      description = "Smoking history (0 = never, 1 = current, 2 = former)",
      units       = "(categorical)",
      type        = "categorical",
      notes       = "Screened (Jin 2025 Resource 1, column TSH) and not retained."
    ),
    ILOC = list(
      description = "Subcutaneous injection location (0 = arm, 1 = stomach, 2 = thigh)",
      units       = "(categorical)",
      type        = "categorical",
      notes       = paste(
        "Screened on absolute subcutaneous bioavailability while diagnosing the AMES under-prediction.",
        "Jin 2025 Discussion: 'Although the injection site covariate reduced the objective function",
        "value, it did not significantly improve the predictive performance for AMES.' Not retained;",
        "the AMES study effect on Fa1 (STUDY_AMES) was selected instead."
      )
    ),
    HEALTHY = list(
      description = "Healthy-volunteer indicator, 1 = healthy subject, 0 = patient with asthma",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened on Fa1, CL and V2 while diagnosing the AMES under-prediction and not retained.",
        "Jin 2025 Discussion: the healthy-status effect, 'likely driven by the larger sample size of",
        "AMES (180 subjects) compared with the D3250C00034 study (36 subjects), improved the",
        "predictive performance for the AMES study but negatively impacted D3250C00034'. The",
        "study-specific STUDY_AMES effect on Fa1 was selected instead."
      )
    ),
    ASSAY = list(
      description = "Bioanalytical assay used (1 = Meso Scale Discovery, 2 = ELISA)",
      units       = "(categorical)",
      type        = "categorical",
      notes       = paste(
        "Screened (Jin 2025 Resource 1, column Assay). Not retained as a parameter covariate; assay /",
        "study-era differences are instead absorbed by the three study-stratified residual-error",
        "magnitudes (see STUDY_PHASE3). The three new studies used Meso Scale Discovery with an LLOQ",
        "of 3.86 ng/mL; one legacy study used an assay with an LLOQ of 60 ng/mL (Jin 2025 Methods 2.3)."
      )
    ),
    SPECIMEN = list(
      description = "Matrix the PK sample was measured in (1 = serum, 2 = plasma)",
      units       = "(categorical)",
      type        = "categorical",
      notes       = "Screened (Jin 2025 Resource 1, column Medium) and not retained."
    ),
    CONMED_MACROLIDE = list(
      description = "Concomitant macrolide use",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened (Jin 2025 Resource 1, column CMACR; Resource 2 legacy covariate list) and not retained."
    ),
    CONMED_MONTELUKAST = list(
      description = "Concomitant montelukast use",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened (Jin 2025 Resource 1, column CMONT; Resource 2 legacy covariate list) and not retained."
    ),
    CONMED_PARACETAMOL = list(
      description = "Concomitant paracetamol (acetaminophen) use",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened (Jin 2025 Resource 1, column CPARA; Resource 2 legacy covariate list) and not retained."
    ),
    CONMED_PPI = list(
      description = "Concomitant proton-pump-inhibitor use",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened (Jin 2025 Resource 1, column CPPI; Resource 2 legacy covariate list) and not retained."
    ),
    CONMED_THEOPHYLLINE = list(
      description = "Concomitant theophylline / aminophylline use",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened on the PK parameters (Jin 2025 Resource 1, column CTHEO) and not retained.",
        "Theophylline co-medication IS a retained covariate on baseline FEV1 in the FEV1",
        "exposure-response layer -- see modellib('Jin_2025_benralizumab_fev1')."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 2855L,
    n_studies      = 12L,
    age_range      = "12-75 years (adults 18-75; adolescents 12-17; no participants aged <12 years were enrolled)",
    age_median     = "48.7 +/- 13 years (mean +/- SD) in adults and 14.4 +/- 1.74 years in adolescents; median not published",
    weight_range   = "40-204 kg",
    weight_median  = "77.5 +/- 18.9 kg (mean +/- SD) in adults and 60.6 +/- 20.6 kg in adolescents; median not published",
    sex_female_pct = 60.5,
    race_ethnicity = c(Asian = 20.7, `Asian (Chinese subset)` = 9.8, `non-Asian` = 79.3),
    disease_state  = paste(
      "Pooled: severe, uncontrolled eosinophilic asthma (phase II/III patients) and healthy",
      "volunteers (phase I). Benralizumab is given as add-on maintenance therapy."
    ),
    dose_range     = paste(
      "0.0003-3 mg/kg intravenously as a single dose (two phase I and one phase II study);",
      "2-200 mg subcutaneously Q4W or Q8W, with the first three doses Q4W (two phase II studies);",
      "10-100 mg subcutaneously as a single dose (two phase I studies); and 30 mg subcutaneously",
      "Q4W or Q8W (five phase III studies)."
    ),
    regions        = paste(
      "Multi-regional. 12 pooled phase I-III studies, including the phase III MIRACLE study",
      "(NCT03186209) in 695 Asian patients aged 12-75 years, a phase I study (NCT03928262) in 36",
      "healthy Han Chinese volunteers, and the phase I AMES study (NCT02968914) in 180 healthy",
      "volunteers. Asian participants were mostly East Asian."
    ),
    notes          = paste(
      "Baseline demographics from Jin 2025 Table 2 (non-Chinese n = 2574, Chinese n = 281, adults",
      "n = 2797, adolescents n = 58). 17,465 serum/plasma benralizumab concentrations from 2855",
      "participants; the nine legacy studies contributed 14,918 observations from 2317 participants.",
      "Chinese participants were 10.0% of the pool and Asian participants 20.7%; the two covariates",
      "are correlated with Phi = 64.2%. Mean body weight was 14.5 kg lower in Chinese (64.1 kg) than",
      "in non-Chinese (78.6 kg) participants. The 58 adolescents came from SIROCCO (29), CALIMA (28)",
      "and MIRACLE (1, from the Philippines); only 2 adolescents were Asian and none were Chinese, so",
      "the paediatric and Chinese-adolescent conclusions rest on simulation rather than observed data.",
      "0.832% of observations were below the limit of quantification and were excluded from the fit.",
      "The model predicts a terminal half-life of 15.5 days (Jin 2025 Results 3.2). Estimation used",
      "NONMEM 7.5.1 with MU referencing.",
      "Paediatric simulations sampled body weight from the Chinese growth references of Zong and Li",
      "(2013) via the R package childsds; children weighing <35 kg received 10 mg Q8W and those",
      ">=35 kg received 30 mg Q8W."
    )
  )

  ini({
    # ------------------------------------------------------------------------
    # Structural disposition parameters. All values are the final-model
    # estimates in Jin 2025 Resource 10 ("Parameter estimates for the final
    # model"), which is the "Global legacy model (full data) + AMES effect Fa1
    # + ASIAN effect CL" column of Resource 11. Reference covariate values:
    # 70 kg body weight, non-Asian, ADA-negative, no study or dose effect.
    # ------------------------------------------------------------------------
    lcl <- log(0.269);  label("Systemic elimination clearance CL (L/day)")                 # Jin 2025 Resource 10, "CL (L/day)" = 0.269 (RSE 1.93%)
    lvc <- log(3.02);   label("Central volume of distribution V2 (L)")                     # Jin 2025 Resource 10, "V2 (L)" = 3.02 (RSE 3.7%)
    lq  <- log(1.05);   label("Intercompartmental clearance Q2 (L/day)")                   # Jin 2025 Resource 10, "Q2 (L/day)" = 1.05 (RSE 4.62%)
    lvp <- log(2.67);   label("Peripheral volume of distribution V3 (L)")                  # Jin 2025 Resource 10, "V3 (L)" = 2.67 (RSE 4.02%)

    # Jin 2025 parameterises subcutaneous absorption by its HALF-LIFE
    # (KAThalf = 3.02 days, RSE 5.22%), not by the rate constant, so the
    # canonical lka is derived: ka = log(2) / KAThalf = 0.2295 /day. Because a
    # lognormal random effect on a half-life is exactly a lognormal random
    # effect of the same magnitude on the corresponding rate constant (with the
    # sign of eta mirrored, which leaves the distribution unchanged), the
    # tabulated omega for KAThalf transfers unaltered to etalka below.
    lka <- log(log(2) / 3.02);  label("First-order subcutaneous absorption rate ka (1/day)")  # Jin 2025 Resource 10, "KAThalf (day)" = 3.02 -> ka = log(2)/3.02

    # ------------------------------------------------------------------------
    # Absolute subcutaneous bioavailability. Jin 2025 estimates three separate
    # absolute bioavailability FRACTIONS (not relative shifts), each with its
    # own IIV: a reference value, a study MI-CP220 value, and an AMES study
    # value. A subject belongs to exactly one of the three strata, so the
    # strata are mutually exclusive in model(). Stratum-suffixed parameter
    # names follow parameter-names.md "Stratum-suffixed parameters".
    # ------------------------------------------------------------------------
    lfdepot         <- log(0.539);  label("Absolute subcutaneous bioavailability Fa1, reference studies (fraction)")  # Jin 2025 Resource 10, "Fa1 (fraction)" = 0.539 (RSE 2.28%)
    lfdepot_micp220 <- log(0.457);  label("Absolute subcutaneous bioavailability, study MI-CP220 (fraction)")         # Jin 2025 Resource 10, "Fa1S220 (fraction)" = 0.457 (RSE 3.69%)
    lfdepot_ames    <- log(0.688);  label("Absolute subcutaneous bioavailability, AMES study (fraction)")             # Jin 2025 Resource 10, "FaslS30 (fraction)" = 0.688 (RSE 2.67%)

    # ------------------------------------------------------------------------
    # Covariate effects. Continuous covariates are multiplicative power terms
    # normalised to a reference value and categorical covariates are
    # exp(beta * indicator), per the covariate-model equations printed in
    # Jin 2025 Resource 2 ("Covariate modeling approach").
    # ------------------------------------------------------------------------
    e_wt_cl    <- 0.849;   label("Power exponent of WT/70 on CL (unitless)")                                   # Jin 2025 Resource 10, "Beta_CL, BWGT (kg)" = 0.849 (RSE 4.03%)
    e_wt_vc    <- 0.799;   label("Power exponent of WT/70 on V2 (unitless)")                                   # Jin 2025 Resource 10, "Beta_V2, BWGT (kg)" = 0.799 (RSE 13.2%)
    e_wt_vp    <- 0.639;   label("Power exponent of WT/70 on V3 (unitless)")                                   # Jin 2025 Resource 10, "Beta_V3, BWGT (kg)" = 0.639 (RSE 15.4%)
    e_asian_cl <- 0.0952;  label("Log-scale effect of Asian race on CL (unitless; +9.99% on CL)")              # Jin 2025 Resource 10, "Beta_CL, ASIAN_1" = 0.0952 (RSE 18.2%), Tvalue +9.99%
    e_ada_cl   <- 0.762;   label("Log-scale effect of ADA positivity on CL (unitless; +114% on CL)")           # Jin 2025 Resource 10, "RCLADA" = 0.762 (RSE 1.7%), Tvalue +114%

    e_dosehigh_fdepot <- -0.554;  label("Log-scale effect of the 200 mg dose on Fa1 (unitless; -42.5% on Fa1)")  # Jin 2025 Resource 10, "Rfa1Dose" = -0.554 (RSE 5.74%), Tvalue -42.5%

    # ------------------------------------------------------------------------
    # Inter-individual variability. Jin 2025 Resource 10 footnote c: "Omega
    # values in this table are presented as standard deviations". nlmixr2's
    # `~` takes a VARIANCE, so each tabulated SD is squared below. Footnote b:
    # "Inter-individual variability parameters were all lognormal, except
    # omega (Rfa1Dose) and omega (RCLADA)" -- those two are therefore encoded
    # as ADDITIVE (normal) random effects on the corresponding log-scale
    # covariate coefficients rather than as lognormal effects on a parameter.
    # ------------------------------------------------------------------------
    etalcl ~ 0.053824              # Jin 2025 Resource 10: omega(CL) SD 0.232 (RSE 3.57%, shrinkage 29.2%) -> 0.232^2; reported 23.5% CV
    etalvc ~ 0.077284              # Jin 2025 Resource 10: omega(V2) SD 0.278 (RSE 7.38%, shrinkage 66.2%) -> 0.278^2; reported 28.3% CV
    etalq  ~ fixed(0.007921)       # Jin 2025 Resource 10: omega(Q2) SD 0.089, held constant by the authors (shrinkage 91.9%) -> 0.089^2; reported 8.9% CV
    etalvp ~ 0.184041              # Jin 2025 Resource 10: omega(V3) SD 0.429 (RSE 4.26%, shrinkage 45.2%) -> 0.429^2; reported 45.0% CV
    etalka ~ 0.499849              # Jin 2025 Resource 10: omega(KAThalf) SD 0.707 (RSE 4.7%, shrinkage 57.1%) -> 0.707^2; reported 80.5% CV. Transfers unchanged from KAThalf to ka (see lka).

    etalfdepot         ~ 0.075076  # Jin 2025 Resource 10: omega(Fa1) SD 0.274 (RSE 4.3%, shrinkage 50.2%) -> 0.274^2; reported 27.9% CV
    etalfdepot_micp220 ~ 0.139129  # Jin 2025 Resource 10: omega(Fa1S220) SD 0.373 (RSE 9.68%, shrinkage 77.9%) -> 0.373^2; reported 38.6% CV
    etalfdepot_ames    ~ 0.020736  # Jin 2025 Resource 10: omega(FalS30) SD 0.144 (RSE 14.4%, shrinkage 84.5%) -> 0.144^2; reported 14.5% CV

    # Additive (normal) random effects on two covariate coefficients, both held
    # fixed by the authors and both essentially degenerate (shrinkage 99.3%).
    etae_dosehigh_fdepot ~ fixed(0.000841)  # Jin 2025 Resource 10: omega(Rfa1Dose) SD 0.029, held constant by the authors (shrinkage 99.3%) -> 0.029^2; normal, not lognormal (footnote b)
    etae_ada_cl          ~ fixed(0.001296)  # omega(RCLADA) SD 0.036, held constant by the authors -> 0.036^2; normal, not lognormal (footnote b). Jin 2025 Resource 8 and Resource 11 both report 0.036 (FIX) in every model column; Resource 10 prints 0.029 for this row, which duplicates the Rfa1Dose row directly above it and is treated as a transcription slip -- see the vignette Errata.

    # ------------------------------------------------------------------------
    # Residual error. Jin 2025 reports three ADDITIVE-ON-THE-LOG-SCALE residual
    # errors in log(ng/mL), one per study stratum (Resource 10 footnote:
    # "ADD1 additive error (all observations except phase III and study
    # MI-CP220), ADD2 additive error (study MI-CP220), ADD3 additive error
    # (SIROCCO, CALIMA, ZONDA, BISE, and MIRACLE studies)"). Additive error on
    # log(concentration) is exactly nlmixr2's lnorm() error model, so the
    # tabulated SDs are used unchanged and are dimensionless on the log scale
    # (they do NOT depend on whether concentration is carried in ng/mL or
    # mg/L). Because nlmixr2 selects a residual error by ENDPOINT rather than
    # by a covariate value, the three strata are encoded as three endpoints
    # over the same underlying prediction; assign each concentration record to
    # the endpoint matching its study (see covariateData$STUDY_PHASE3).
    # ------------------------------------------------------------------------
    expSd         <- 0.367;  label("Log-scale additive residual error, phase III studies (log(ng/mL))")                 # Jin 2025 Resource 10, "Error_ADD3, log(ng/mL)" = 0.367 (RSE 0.737%) -- SIROCCO, CALIMA, ZONDA, BISE, MIRACLE
    expSd_Cc_early   <- 0.175;  label("Log-scale additive residual error, early phase I/II studies (log(ng/mL))")          # Jin 2025 Resource 10, "Error_ADD1, log(ng/mL)" = 0.175 (RSE 1.64%) -- all observations except phase III and MI-CP220
    expSd_Cc_micp220 <- 0.549;  label("Log-scale additive residual error, study MI-CP220 (log(ng/mL))")                    # Jin 2025 Resource 10, "Error_ADD2, log(ng/mL)" = 0.549 (RSE 2.17%) -- study MI-CP220
  })
  model({
    # --- 1. Individual disposition parameters -------------------------------
    # Covariate structure per Jin 2025 Resource 2: continuous covariates enter
    # as (COV/REF)^beta and categorical covariates as exp(beta * indicator),
    # multiplied together. Reference subject: 70 kg, non-Asian, ADA-negative.
    cl <- exp(lcl + etalcl) *
          (WT / 70)^e_wt_cl *
          exp(e_asian_cl * RACE_ASIAN) *
          exp((e_ada_cl + etae_ada_cl) * ADA_POS)
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vp
    q  <- exp(lq  + etalq)
    ka <- exp(lka + etalka)

    # --- 2. Absolute subcutaneous bioavailability --------------------------
    # Exactly one of the three study strata applies to a given subject. The
    # 200 mg dose effect then multiplies whichever stratum value was selected.
    # Each stratum value is built on its own simple line so that every eta stays
    # mu-referenced (rxode2 warns and falls back to non-mu referencing if a
    # mu-referenced expression is buried inside a compound expression).
    fdepotref     <- exp(lfdepot + etalfdepot)
    fdepotmicp220 <- exp(lfdepot_micp220 + etalfdepot_micp220)
    fdepotames    <- exp(lfdepot_ames + etalfdepot_ames)
    fdepotstudy   <- fdepotref * (1 - STUDY_MICP220) * (1 - STUDY_AMES) +
                     fdepotmicp220 * STUDY_MICP220 +
                     fdepotames * STUDY_AMES
    fdepot <- fdepotstudy * exp((e_dosehigh_fdepot + etae_dosehigh_fdepot) * DOSE_HIGH)

    # --- 3. Micro-constants ------------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # --- 4. Two-compartment ODE system with a first-order SC depot ----------
    # Intravenous doses in the phase I/II cohorts are given directly into
    # `central`, which bypasses `depot` and therefore also bypasses fdepot --
    # this is what makes Fa1 an ABSOLUTE (rather than relative) subcutaneous
    # bioavailability, as Jin 2025 Resource 10 labels it.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # --- 5. Bioavailability -------------------------------------------------
    f(depot) <- fdepot

    # --- 6. Observation and study-stratified residual error -----------------
    Cc <- central / vc

    # Cc_early and Cc_micp220 are the SAME predicted serum concentration as Cc;
    # they exist only so that each study stratum can carry its own residual
    # error magnitude, which nlmixr2 keys off the endpoint. Simulation users
    # who do not need the stratified error can read the Cc column alone.
    Cc_early   <- Cc
    Cc_micp220 <- Cc

    Cc         ~ lnorm(expSd)
    Cc_early   ~ lnorm(expSd_Cc_early)
    Cc_micp220 ~ lnorm(expSd_Cc_micp220)
  })
}
