Yan_2019_benralizumab <- function() {
  description <- "Two-compartment population PK model of benralizumab (anti-IL-5R alpha) with first-order subcutaneous absorption and first-order elimination in adult and adolescent patients with asthma (Yan 2019), pooling nine phase I-III studies; body weight on CL/Vc/Vp, anti-drug antibody on CL, a study-specific absolute subcutaneous bioavailability for the phase IIb study MI-CP220, and study-stratified log-scale residual error"
  reference <- paste(
    "Yan L, Wang B, Chia YL, Roskos LK.",
    "Population pharmacokinetic modeling of benralizumab in adult and adolescent",
    "patients with asthma.",
    "Clin Pharmacokinet. 2019;58(7):943-958. doi:10.1007/s40262-019-00738-4.",
    "Parameter values are the 'Final updated model (Model 9)' column of Table 5,",
    "the model refitted to all nine studies after ZONDA was appended to the",
    "analysis dataset. Superseded for the Asian and paediatric populations by",
    "Jin Y, Guiastrennec B, Stuke M, et al. Clin Pharmacokinet. 2025;64:1233-1245;",
    "doi:10.1007/s40262-025-01538-9 -- see modellib('Jin_2025_benralizumab').",
    sep = " "
  )
  vignette <- "Yan_2019_benralizumab"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(analyte = "benralizumab", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "benralizumab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "benralizumab", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power (multiplicative) effects normalised to 70 kg on the three disposition parameters:",
        "(WT/70)^0.807 on CL, (WT/70)^0.803 on Vc and (WT/70)^0.528 on Vp (Yan 2019 Table 5,",
        "'Final updated model (Model 9)' column). The 70 kg reference is stated twice: Yan 2019",
        "Methods 2.5.6 -- 'the covariate X normalized by a reference value Xref (70 kg for body",
        "weight and approximate median for other covariates)' -- and the Figure 4 caption, which",
        "prints the prediction curve as 'CL.(body weight/70)^0.807'. Yan 2019 Discussion calls the",
        "CL effect 'nearly allometric'. A body-weight power on Q was also carried in the full model",
        "but is reported as 0 with no confidence interval in Table 5, i.e. the WAM procedure of",
        "Methods 2.5.7 dropped it; Q therefore has NO weight effect in the final model and no",
        "zero-valued exponent is encoded here. Baseline value, time-fixed per subject. Yan 2019",
        "Table 2 gives 79.14 +/- 19.74 kg (mean +/- SD), median 77 kg, range 40-204.4 kg over the",
        "pooled 2316 patients with a recorded weight. Source column WT."
      ),
      source_name = "WT"
    ),
    ADA_POS = list(
      description = "Anti-drug antibody status, 1 = positive at any assessed visit, 0 = negative",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (anti-drug antibody negative)",
      notes = paste(
        "Multiplicative effect on clearance: CL is 2.24-fold higher in ADA-positive patients",
        "(Yan 2019 Table 5 row 'ADAs on CL, fraction' = 2.24, 90% CI 2.18-2.30). The abstract and",
        "Discussion state the same effect as a percentage -- 'the presence of ADAs increased",
        "benralizumab CL by 124%' -- which fixes 2.24 as the back-transformed exp(theta) factor of",
        "the categorical fractional-change model in Methods 2.5.6 Eq. 4, not as a fitted log-scale",
        "coefficient. It is encoded below as e_ada_cl = log(2.24). Adding ADA status to the model",
        "gave a >4000-unit drop in objective function (Yan 2019 Results 3.3 and Table 4, Model 7),",
        "making it the single most influential covariate in the analysis. ADA status was assessed at",
        "each predesignated sampling visit and is therefore time-varying in the source dataset;",
        "362 of 2317 patients (15.62%) were ever positive (Yan 2019 Table 3). Yan 2019 uses a plain",
        "positive/negative immunogenicity status -- the sibling modellib('Wang_2017_benralizumab')",
        "instead uses a HIGH-TITER (>=400) definition and the two are not interchangeable;",
        "modellib('Jin_2025_benralizumab') uses the same plain positive/negative definition as here.",
        "Source column ADA."
      ),
      source_name = "ADA"
    ),
    STUDY_MICP220 = list(
      description = "Study MI-CP220 indicator, 1 = patient enrolled in the phase IIb benralizumab study MI-CP220 (NCT01238861), 0 = any of the other eight pooled studies",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (any of the other eight pooled studies)",
      notes = paste(
        "Selects a study-specific absolute subcutaneous bioavailability of 0.490 in place of the",
        "reference 0.589 (Yan 2019 Table 5 rows 'F' and 'Change in F with study CP220, fraction';",
        "'CP220' is the paper's short form of MI-CP220 / NCT01238861). Yan 2019 Results 3.2:",
        "'The estimated F1 and BPV of F1 for study NCT01238861 were different from the other",
        "studies.' The 0.490 is the study's ABSOLUTE bioavailability fraction, not a multiplier on",
        "0.589 -- three independent checks agree: (i) the same row in the base-model column gives",
        "0.569 against a reference F of 0.671, a ratio of 0.848, and the successor analysis",
        "modellib('Jin_2025_benralizumab') fits the same MI-CP220 data and reports Fa1S220 / Fa1 =",
        "0.457 / 0.539 = 0.848, an exact match that a 0.569 multiplier could not produce;",
        "(ii) the final-model ratio 0.490 / 0.589 = 0.832 stays in the same place, whereas the",
        "multiplier reading would move the relative bioavailability from 0.569 to 0.490 between two",
        "fits of overlapping data; (iii) the stratum carries its own separately tabulated IIV",
        "('etaF (study CP220), %CV' = 35.0 versus 17.1 for the reference F), which is the signature",
        "of a separately estimated absolute fraction rather than of a shift applied to a shared",
        "parameter. MI-CP220 is also the single study with its own residual-error magnitude",
        "(54.5%CV versus 25.0%CV for the other early studies), encoded here as the separate",
        "Cc_micp220 endpoint. Time-fixed per subject. Derived from the study identifier column."
      ),
      source_name = "STUDY"
    )
  )

  # Covariates that Yan 2019 screened but did not retain in the final model.
  # Documented for provenance only; deliberately absent from model().
  covariatesDataExcluded <- list(
    STUDY_PHASE3 = list(
      description = "Phase III study indicator, 1 = patient enrolled in SIROCCO, CALIMA, ZONDA or BISE, 0 = patient enrolled in one of the five phase I/II studies",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "DATASET-CONSTRUCTION GUIDANCE, not a model covariate. It records which of the three",
        "study-stratified residual-error magnitudes applies to a concentration record, and therefore",
        "which of the three model endpoints the record belongs to. Yan 2019 Results 3.2: 'the",
        "variance of residual error (sigma^2) was different for early stage studies (phase I-IIa,",
        "studies NCT00512486, NCT00659659, NCT00768079 and NCT00783289), phase IIb (NCT01238861),",
        "and phase III studies (SIROCCO, CALIMA, and BISE)'; ZONDA joined the phase III stratum when",
        "it was appended for the final updated model, as the Table 5 row label",
        "'Proportional error (studies SIROCCO, CALIMA, ZONDA, BISE)' shows. Records with",
        "STUDY_PHASE3 = 1 are the Cc endpoint (expSd = 0.367); records with STUDY_MICP220 = 1 are",
        "the Cc_micp220 endpoint (expSd_Cc_micp220 = 0.545); all remaining records are the Cc_early",
        "endpoint (expSd_Cc_early = 0.250). It is deliberately NOT referenced in model(), because",
        "nlmixr2 selects a residual error by endpoint (cmt/dvid) rather than by a covariate value;",
        "the stratum is expressed by which endpoint a row is assigned to."
      )
    ),
    AGE = list(
      description = "Baseline subject age",
      units = "years",
      type = "continuous",
      notes = paste(
        "Screened both as a continuous covariate and as an adult/adolescent age-group flag",
        "(Yan 2019 Results 3.1) and not retained. Yan 2019 Discussion: 'Benralizumab mean exposure",
        "decreased slightly from adolescents (12-17 years of age) to older adults (65-75 years of",
        "age). The terminal half-life of benralizumab increased slightly with age, with median",
        "estimates of approximately 14 days in adolescents, 15 days in adults (18-64 years of age),",
        "and 17 days in older adults.' Body weight carries the adolescent size effect instead.",
        "Yan 2019 Table 2 gives 48.48 +/- 13.88 years, median 50 (12-75); only 57 of 2317 patients",
        "(2.46%) were adolescents, all from SIROCCO and CALIMA (Table 3)."
      )
    ),
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Screened on CL and Vc and not retained. Yan 2019 Results 3.3: 'Using the WAM algorithm,",
        "the effect of sex, race-Asian, and race-black on CL was <10%... These covariate effects",
        "were also considered nonmeaningful for Vc.' 1482 of 2317 patients (63.96%) were female",
        "(Yan 2019 Table 3)."
      )
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator, 1 = Asian, 0 = any other race category",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Screened on CL and Vc and not retained: the WAM effect on CL was <10% (Yan 2019",
        "Results 3.3). 227 of 2317 patients (9.8%) were Asian (Table 3). The successor analysis",
        "modellib('Jin_2025_benralizumab'), which added a phase III study in Asian patients, DID",
        "retain an Asian effect on CL (+9.99%) even though it likewise falls below the 10%",
        "clinical-relevance threshold, because it improved predictive performance for that study."
      )
    ),
    RACE_BLACK = list(
      description = "Black / African American race indicator, 1 = African American, 0 = any other race category",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Screened on CL and Vc and not retained: the WAM effect on CL was <10% (Yan 2019",
        "Results 3.3). 115 of 2317 patients (4.96%) were African American (Table 3)."
      )
    ),
    RACE_OTHER = list(
      description = "Race-category 'Other' indicator, 1 = race reported as Other, 0 = White, African American or Asian",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Screened on CL and not retained. This is the ONE race category whose WAM point estimate",
        "exceeded the 10% threshold: Yan 2019 Results 3.3 reports 'The effect of race-other on CL",
        "was 19%, but the 90% confidence interval precluded a meaningful effect.' 174 of 2317",
        "patients (7.51%) were in this category, 76 of them from the phase IIb study NCT01238861",
        "(Table 3). The reference category was White (1801 patients, 77.73%)."
      )
    ),
    SMOKE = list(
      description = "Smoking history (never / former / current smoker)",
      units = "(categorical)",
      type = "categorical",
      notes = paste(
        "Screened (Yan 2019 Results 3.1 and Table 3) and not retained. Only 26 of 2317 patients",
        "(1.12%) were current smokers -- the phase III protocols excluded current smokers -- so the",
        "current-smoker stratum carried almost no information."
      )
    ),
    ALB = list(
      description = "Baseline serum albumin",
      units = "g/L",
      type = "continuous",
      notes = paste(
        "Screened (Yan 2019 Results 3.1) and not retained. Recorded for only 497 of 2317 patients",
        "(Yan 2019 Table 2: 44.09 +/- 2.89 g/L, median 44, range 35-52), because the two largest",
        "phase III studies contributed 25 albumin values between them."
      )
    ),
    CRCL = list(
      description = "Baseline creatinine clearance (Cockcroft-Gault)",
      units = "mL/min",
      type = "continuous",
      notes = paste(
        "The only covariate with a visible eta trend, and deliberately still excluded. Yan 2019",
        "Results 3.3: 'The effect of CRCL on CL was statistically significant; an estimate of 0 was",
        "not precluded. There was a general trend of increasing effect on CL with greater CRCL.",
        "However, the apparent correlation of CRCL and benralizumab CL could be an artifact as both",
        "were influenced by body weight (heavier patients tended to have greater CRCL and CL). A",
        "posterior ad hoc analysis confirmed that renal function (estimated glomerular filtration",
        "rate [eGFR]) did not influence benralizumab CL.' Figure 5a shows benralizumab CL rising",
        "from 0.23 to 0.33 L/day across CRCL strata, while Figure 5b shows it flat at 0.30-0.31",
        "L/day across eGFR strata -- eGFR being the body-size-normalised marker. Yan 2019 Table 2:",
        "111.46 +/- 36.15 mL/min, median 105.71 (2.41-349.17). Note that the canonical CRCL column",
        "is BSA-normalised (mL/min/1.73 m^2); Yan 2019 reports raw Cockcroft-Gault mL/min, which is",
        "precisely why the effect is confounded with body weight."
      )
    ),
    EGFR = list(
      description = "Baseline estimated glomerular filtration rate",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      notes = paste(
        "Evaluated post hoc specifically to break the CRCL/body-weight confounding described in the",
        "CRCL entry above, and not retained -- Yan 2019 Discussion: 'A post hoc analysis using eGFR",
        "as a marker demonstrated that renal function had no effect on the PK of benralizumab.'",
        "Expected for a 150 kDa IgG1 monoclonal antibody, which is not renally eliminated.",
        "Yan 2019 Table 2: 97.49 +/- 33.06 mL/min/1.73 m^2, median 97.78 (0.87-251.62)."
      )
    ),
    ALP = list(
      description = "Baseline alkaline phosphatase",
      units = "ukat/L",
      type = "continuous",
      notes = "Hepatic marker screened (Yan 2019 Results 3.1 and Table 2: 1.24 +/- 0.52 ukat/L) and not retained."
    ),
    ALT = list(
      description = "Baseline alanine aminotransferase",
      units = "ukat/L",
      type = "continuous",
      notes = "Hepatic marker screened (Yan 2019 Results 3.1 and Table 2: 0.39 +/- 0.25 ukat/L) and not retained."
    ),
    AST = list(
      description = "Baseline aspartate aminotransferase",
      units = "ukat/L",
      type = "continuous",
      notes = "Hepatic marker screened (Yan 2019 Results 3.1 and Table 2: 0.35 +/- 0.15 ukat/L) and not retained."
    ),
    TBIL = list(
      description = "Baseline total bilirubin",
      units = "umol/L",
      type = "continuous",
      notes = paste(
        "Hepatic marker screened (Yan 2019 Results 3.1 and Table 2: 8 +/- 4.2 umol/L) and not",
        "retained. Yan 2019 Discussion: 'IgG mAbs are not primarily cleared via the hepatic pathway;",
        "therefore, change in hepatic function is not expected to influence benralizumab CL.'"
      )
    ),
    EOS = list(
      description = "Baseline blood eosinophil count",
      units = "cells/uL",
      type = "continuous",
      notes = paste(
        "Screened on the PK parameters and not retained. This is a mechanistically load-bearing",
        "negative: benralizumab is an IL-5R alpha-directed cytolytic antibody, so a target-mediated",
        "disposition pathway would show up as a baseline-eosinophil effect on CL. Yan 2019",
        "Discussion: 'The lack of any trend in CL across groups suggests that eosinophil count does",
        "not have a meaningful impact on benralizumab CL, as evidenced by the absence of a nonlinear",
        "IL-5R-mediated elimination pathway of benralizumab in humans', which Results 3.2 attributes",
        "to 'the rapid depletion of circulating IL-5R-expressing eosinophils following benralizumab",
        "administration'. Figure 5c shows CL flat at 0.30-0.31 L/day across four baseline-eosinophil",
        "strata (<150, 150-300, 300-449, >=450 cells/uL). This is why the model below is linear."
      )
    ),
    ILOC = list(
      description = "Subcutaneous injection location (arm, stomach or thigh)",
      units = "(categorical)",
      type = "categorical",
      notes = paste(
        "Evaluated ad hoc on bioavailability rather than through the eta plots (Yan 2019",
        "Methods 2.5.6: 'The effects of injection site (arm, stomach, or thigh) and antidrug",
        "antibody (ADA) status could not be evaluated by covariate eta plots and were evaluated in",
        "an ad hoc fashion') and not retained. Yan 2019 Discussion: 'The estimated subcutaneous",
        "bioavailability of benralizumab administered to the stomach or thigh was approximately 9%",
        "greater than the bioavailability after administration to the upper arm; however, the",
        "magnitude of the effect was not considered clinically relevant.' ADA status, evaluated by",
        "the same ad hoc route, WAS retained -- see covariateData$ADA_POS."
      )
    ),
    CONMED_MONTELUKAST = list(
      description = "Concomitant montelukast use",
      units = "(binary)",
      type = "binary",
      notes = "Screened (Yan 2019 Results 3.1 and Table 3: 601 of 2317 patients, 25.94%) and not retained."
    ),
    CONMED_PARACETAMOL = list(
      description = "Concomitant paracetamol (acetaminophen) use",
      units = "(binary)",
      type = "binary",
      notes = "Screened (Yan 2019 Results 3.1 and Table 3: 106 of 2317 patients, 4.57%) and not retained."
    ),
    CONMED_PPI = list(
      description = "Concomitant proton-pump-inhibitor use",
      units = "(binary)",
      type = "binary",
      notes = "Screened (Yan 2019 Results 3.1 and Table 3: 382 of 2317 patients, 16.49%) and not retained."
    ),
    CONMED_MACROLIDE = list(
      description = "Concomitant macrolide use",
      units = "(binary)",
      type = "binary",
      notes = "Screened (Yan 2019 Results 3.1 and Table 3: 12 of 2317 patients, 0.52%) and not retained."
    ),
    CONMED_THEOPHYLLINE = list(
      description = "Concomitant theophylline / aminophylline use",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Screened (Yan 2019 Results 3.1 and Table 3) and not retained. Yan 2019 Discussion groups it",
        "with the other four: the commonly used small-molecule co-medications 'had no effect on",
        "benralizumab CL'."
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 2317L,
    n_studies = 9L,
    age_range = "12-75 years (adults 18-75, n = 2260, 97.54%; adolescents 12-17, n = 57, 2.46%)",
    age_median = "50 years (range 12-75); mean 48.48 +/- 13.88 years",
    weight_range = "40-204.4 kg",
    weight_median = "77 kg (range 40-204.4); mean 79.14 +/- 19.74 kg over the 2316 patients with a recorded weight",
    sex_female_pct = 63.96,
    race_ethnicity = c(
      White = 77.73,
      `African American` = 4.96,
      Asian = 9.8,
      Other = 7.51
    ),
    disease_state = paste(
      "Asthma. The five phase I/II studies enrolled patients with mild atopic, eosinophilic or",
      "uncontrolled asthma; the four phase III studies (SIROCCO, CALIMA, ZONDA, BISE) enrolled",
      "patients with uncontrolled or mild-to-moderate persistent asthma on inhaled",
      "corticosteroids/long-acting beta2-agonists, ZONDA additionally requiring long-term oral",
      "corticosteroid therapy. No healthy volunteers. Benralizumab is given as add-on maintenance",
      "therapy for the eosinophilic phenotype."
    ),
    dose_range = paste(
      "Per Yan 2019 Table 1. Intravenous, single dose: 0.0003-3 mg/kg (phase I NCT00512486),",
      "1.0 mg/kg (phase I NCT00659659) and 0.3 or 1.0 mg/kg (phase II NCT00768079).",
      "Subcutaneous: 100 and 200 mg Q4W for 8 weeks (NCT00659659); 25, 100 or 200 mg Q4W for",
      "8 weeks (phase II NCT00783289); 2, 20 or 100 mg Q4W for 8 weeks then Q8W for 32 weeks",
      "(phase IIb NCT01238861); and 30 mg either Q4W throughout or Q4W for the first three doses",
      "then Q8W (the four phase III studies), plus three doses at weeks 0, 4 and 8 in BISE.",
      "Methods 2.3 summarises the range as 0.0003-3 mg/kg intravenously and 2-200 mg",
      "subcutaneously. 80% of patients received 30 mg subcutaneously and 11% received 100 mg;",
      "the 2 mg subcutaneous arm of NCT01238861 was excluded from the analysis."
    ),
    regions = "Multi-regional; nine AstraZeneca/MedImmune trials, the largest being SIROCCO (NCT01928771) and CALIMA (NCT01914757).",
    immunogenicity = "362 of 2317 patients (15.62%) were anti-drug antibody positive at some assessed visit.",
    notes = paste(
      "14,918 quantifiable serum (plasma for NCT00659659) benralizumab concentrations from 2317",
      "patients across all nine studies support the final updated model (Yan 2019 Results 3.5). The",
      "preceding initial model, from which the base-model column of Table 5 comes, used 14,106",
      "concentrations from 2174 patients across eight studies; ZONDA was held out for external",
      "validation (832 concentrations from 143 patients, Yan 2019 Results 3.1 and 3.4) and then",
      "appended. Excluded before fitting: all placebo patients (n = 1227) and all pre-first-dose",
      "samples; the 589 samples from the 81 patients on benralizumab 2 mg subcutaneously in",
      "NCT01238861, because of low and variable concentrations; and 75 concentrations with",
      "unrealistic values or unclear dosing records. Fewer than 10% of records were below the limit",
      "of quantification, so no BLQ imputation method was applied. Assay LLOQ was 60 ng/mL for",
      "NCT00512486 and 3.86 ng/mL elsewhere. Estimation used NONMEM 7.3 with MU referencing and",
      "expectation-maximization methods (ITS/SAEM with importance sampling) rather than the planned",
      "FOCEI, which 'occasionally led to convergence difficulties'; rerunning the pre-validation",
      "final model under FOCEI changed the objective function by 0.45 units. Covariate selection",
      "used a full model followed by Wald's Approximation Method with a 10.83 penalty (p = 0.001).",
      "The paper reports a terminal elimination half-life of approximately 15.5 days.",
      "Yan 2019 Tables 2 and 3 print 'NCT00768079' twice in the column headers; the fourth column",
      "(n = 19) is NCT00783289, the phase II study named in Table 1 and in Results 3.2."
    )
  )

  ini({
    # ------------------------------------------------------------------------
    # Structural disposition parameters. Every value is the "Final updated
    # model (Model 9)" column of Yan 2019 Table 5 -- the model refitted to all
    # nine studies after ZONDA was appended (Results 3.5). The base-model
    # column of the same table belongs to Model 5, an earlier eight-study fit
    # with no covariates, and is NOT used here.
    #
    # Reference subject: 70 kg body weight, ADA-negative, enrolled in a study
    # other than MI-CP220.
    # ------------------------------------------------------------------------
    lcl <- log(0.291)
    label("Systemic clearance CL (L/day)") # Yan 2019 Table 5, 'CL, L/day' = 0.291 (90% CI 0.28-0.302)
    lvc <- log(3.13)
    label("Central volume of distribution Vc (L)") # Yan 2019 Table 5, 'Vc, L' = 3.13 (90% CI 2.97-3.31)
    lq <- log(0.738)
    label("Intercompartmental clearance Q (L/day)") # Yan 2019 Table 5, 'Q, L/day' = 0.738 (90% CI 0.679-0.803)
    lvp <- log(2.52)
    label("Peripheral volume of distribution Vp (L)") # Yan 2019 Table 5, 'Vp, L' = 2.52 (90% CI 2.34-2.71)

    # Yan 2019 parameterises subcutaneous absorption by its HALF-LIFE, not by
    # the rate constant: Table 5 row 'ka, half-life; days' = 3.54 (90% CI
    # 3.15-3.99), which Results 3.5 and the Discussion restate as 'an estimated
    # absorption half-life of 3.5 days'. The canonical lka is therefore derived,
    # ka = log(2) / 3.54 = 0.1958 /day. A lognormal random effect on a half-life
    # is exactly a lognormal random effect of the same magnitude on the
    # corresponding rate constant (the sign of eta mirrors, which leaves a
    # zero-mean normal distribution unchanged), so the tabulated omega for the
    # half-life transfers unaltered to etalka below.
    lka <- log(log(2) / 3.54)
    label("First-order subcutaneous absorption rate ka (1/day)") # Yan 2019 Table 5, 'ka, half-life; days' = 3.54 -> ka = log(2)/3.54

    # ------------------------------------------------------------------------
    # Absolute subcutaneous bioavailability. Yan 2019 estimates two separate
    # absolute bioavailability FRACTIONS, each with its own IIV: a reference
    # value (Table 5 'F' = 0.589, 90% CI 0.565-0.614) and a value for the
    # phase IIb study MI-CP220 (Table 5 'Change in F with study CP220,
    # fraction' = 0.490, 90% CI 0.461-0.522). See covariateData$STUDY_MICP220
    # for why 0.490 is read as the study's absolute fraction rather than as a
    # multiplier on 0.589. Stratum-suffixed parameter names follow
    # parameter-names.md 'Stratum-suffixed parameters' and match the sibling
    # modellib('Jin_2025_benralizumab'), which uses the same _micp220 suffix
    # for the same study.
    #
    # Intravenous doses go straight into `central` and so bypass both `depot`
    # and fdepot, which is what makes F an ABSOLUTE rather than a relative
    # bioavailability. Results 3.2 records that an exponential parameterisation
    # of F fitted slightly better than a logit transform.
    # ------------------------------------------------------------------------
    lfdepot <- log(0.589)
    label("Absolute subcutaneous bioavailability, reference studies (fraction)") # Yan 2019 Table 5, 'F' = 0.589 (90% CI 0.565-0.614)
    lfdepot_micp220 <- log(0.490)
    label("Absolute subcutaneous bioavailability, study MI-CP220 (fraction)") # Yan 2019 Table 5, 'Change in F with study CP220, fraction' = 0.490 (90% CI 0.461-0.522)

    # ------------------------------------------------------------------------
    # Covariate effects. Continuous covariates enter as (COV/REF)^beta and
    # categorical covariates as exp(beta * indicator), per Yan 2019 Eqs. 3 and
    # 4 (Methods 2.5.6). Table 5 prints the categorical effect already
    # back-transformed as a fraction, so e_ada_cl is written as log(2.24).
    # ------------------------------------------------------------------------
    e_wt_cl <- 0.807
    label("Power exponent of WT/70 on CL (unitless)") # Yan 2019 Table 5, 'Body weight on CL, power' = 0.807 (90% CI 0.751-0.864)
    e_wt_vc <- 0.803
    label("Power exponent of WT/70 on Vc (unitless)") # Yan 2019 Table 5, 'Body weight on Vc, power' = 0.803 (90% CI 0.627-0.979)
    e_wt_vp <- 0.528
    label("Power exponent of WT/70 on Vp (unitless)") # Yan 2019 Table 5, 'Body weight on Vp, power' = 0.528 (90% CI 0.351-0.706)
    e_ada_cl <- log(2.24)
    label("Log-scale effect of ADA positivity on CL (unitless; 2.24-fold, i.e. +124%, on CL)") # Yan 2019 Table 5, 'ADAs on CL, fraction' = 2.24 (90% CI 2.18-2.3); abstract '+124%'

    # ------------------------------------------------------------------------
    # Inter-individual variability. Yan 2019 Table 5 reports every IIV as a
    # percent CV. That percentage is 100 * the omega STANDARD DEVIATION, not
    # the exact lognormal CV: Table 4 Model 8 states 'Fixed IIV(Q) = 0.008' --
    # a NONMEM omega VARIANCE -- and Table 5 prints etaQ as 8.94 %CV, which is
    # 100 * sqrt(0.008) = 8.944 exactly, where the exact lognormal CV
    # 100 * sqrt(exp(0.008) - 1) would have printed as 8.96. nlmixr2's `~`
    # takes a VARIANCE, so each tabulated percentage is squared below.
    # ------------------------------------------------------------------------
    etalcl ~ 0.058564 # Yan 2019 Table 5: etaCL 24.2 %CV (90% CI 22.6-25.6) -> 0.242^2
    etalvc ~ 0.059536 # Yan 2019 Table 5: etaVc 24.4 %CV (90% CI 20.7-27.6) -> 0.244^2
    etalq ~ fixed(0.008) # Yan 2019 Table 4 Model 8 holds IIV(Q) at 0.008; Table 5 prints the matching 8.94 %CV with no confidence interval
    etalvp ~ 0.199809 # Yan 2019 Table 5: etaVp 44.7 %CV (90% CI 40.5-48.5) -> 0.447^2
    etalka ~ 0.690561 # Yan 2019 Table 5: etaka (half-life) 83.1 %CV (90% CI 75.5-90.1) -> 0.831^2; transfers unchanged from the half-life to ka (see lka)

    etalfdepot ~ 0.029241 # Yan 2019 Table 5: etaF 17.1 %CV (90% CI 13.5-20) -> 0.171^2
    etalfdepot_micp220 ~ 0.122500 # Yan 2019 Table 5: etaF (study CP220) 35.0 %CV (90% CI 28.8-40.2) -> 0.350^2

    # ------------------------------------------------------------------------
    # Residual error. Yan 2019 Eq. 2 sets out a combined proportional-plus-
    # additive error, but the analysis was run on LOG-TRANSFORMED
    # concentrations -- Results 3.2: 'the PK concentration data were
    # log-transformed to facilitate model development. A log-normal residual
    # error structure was used to increase the influence of lesser
    # concentrations on the estimates' -- and Table 5 reports only the
    # proportional component, in three study strata. Additive error on
    # log(concentration) is exactly nlmixr2's lnorm() error model, so the
    # tabulated percentages are carried over as log-scale standard deviations
    # (dimensionless; they do not depend on whether concentration is held in
    # mg/L or ng/mL). The successor analysis corroborates the reading: the
    # phase III stratum that Yan 2019 prints as 36.7 %CV is tabulated by
    # modellib('Jin_2025_benralizumab') as 'Error_ADD3, log(ng/mL)' = 0.367.
    #
    # Because nlmixr2 selects a residual error by ENDPOINT rather than by a
    # covariate value, the three strata are encoded as three endpoints over the
    # same underlying prediction; assign each concentration record to the
    # endpoint matching its study (see covariatesDataExcluded$STUDY_PHASE3).
    # ------------------------------------------------------------------------
    expSd <- 0.367
    label("Log-scale residual error, phase III studies (unitless)") # Yan 2019 Table 5, 'Proportional error (studies SIROCCO, CALIMA, ZONDA, BISE), %CV' = 36.7 (90% CI 36.2-37.1)
    expSd_Cc_early <- 0.250
    label("Log-scale residual error, phase I-IIa studies (unitless)") # Yan 2019 Table 5, 'Proportional error, %CV' = 25.0 (90% CI 23.8-26.1); Results 3.2 assigns this row to NCT00512486, NCT00659659, NCT00768079 and NCT00783289
    expSd_Cc_micp220 <- 0.545
    label("Log-scale residual error, study MI-CP220 (unitless)") # Yan 2019 Table 5, 'Proportional error (study CP220), %CV' = 54.5 (90% CI 52.5-56.3)
  })
  model({
    # --- 1. Individual disposition parameters -------------------------------
    # Yan 2019 Eq. 3 (continuous, power) and Eq. 4 (categorical, fractional
    # change), both on the log-based parameterisation of Eq. 1. Reference
    # subject: 70 kg, ADA-negative. Body weight has no effect on Q -- the full
    # model carried one, but Table 5 reports the final exponent as 0 with no
    # confidence interval.
    cl <- exp(lcl + etalcl) *
      (WT / 70)^e_wt_cl *
      exp(e_ada_cl * ADA_POS)
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vp
    q <- exp(lq + etalq)
    ka <- exp(lka + etalka)

    # --- 2. Absolute subcutaneous bioavailability --------------------------
    # A patient belongs to exactly one of the two study strata. Each stratum
    # value is built on its own simple line so that every eta stays
    # mu-referenced (rxode2 warns and falls back to non-mu referencing if a
    # mu-referenced expression is buried inside a compound expression).
    fdepotref <- exp(lfdepot + etalfdepot)
    fdepotmicp220 <- exp(lfdepot_micp220 + etalfdepot_micp220)
    fdepot <- fdepotref * (1 - STUDY_MICP220) + fdepotmicp220 * STUDY_MICP220

    # --- 3. Micro-constants ------------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # --- 4. Two-compartment ODE system with a first-order SC depot ----------
    # Yan 2019 Fig. 1. Intravenous doses in the phase I/II cohorts go directly
    # into `central`, bypassing `depot` and therefore also bypassing fdepot.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # --- 5. Bioavailability -------------------------------------------------
    f(depot) <- fdepot

    # --- 6. Observation and study-stratified residual error -----------------
    Cc <- central / vc

    # Cc_early and Cc_micp220 are the SAME predicted serum concentration as Cc;
    # they exist only so that each study stratum can carry its own residual
    # error magnitude, which nlmixr2 keys off the endpoint. Simulation users
    # who do not need the stratified error can read the Cc column alone.
    Cc_early <- Cc
    Cc_micp220 <- Cc

    Cc ~ lnorm(expSd)
    Cc_early ~ lnorm(expSd_Cc_early)
    Cc_micp220 ~ lnorm(expSd_Cc_micp220)
  })
}
