Miano_2024_tacrolimus <- function() {
  description <- "One-compartment population PK model for oral / sublingual tacrolimus during the first 14 postoperative days after lung transplantation (Miano 2024), with first-order absorption at a fixed ka of 4.5 1/h, a power-of-postoperative-day effect on apparent clearance, a negative power-of-hematocrit effect on CL/F, a time-varying bilateral-versus-single lung graft effect on CL/F (separate coefficients for postoperative days 1-3 and 4-14), multiplicative CYP3A5-expresser, voriconazole and amiodarone effects on CL/F, power-of-weight scaling on both CL/F and Vd/F, exponential inter-individual variability on CL/F and Vd/F, and combined additive plus proportional residual error."
  reference <- "Miano TA, Zuppa AF, Feng R, Griffiths S, Kalman L, Oyster M, Cantu E, Yang W, Diamond JM, Christie JD, Scheetz MH, Shashaty MGS. Development and validation of a population pharmacokinetic model to guide perioperative tacrolimus dosing after lung transplantation. JHLT Open. 2024;6:100134. doi:10.1016/j.jhlto.2024.100134"
  vignette <- "Miano_2024_tacrolimus"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Tacrolimus concentrations in Miano 2024 were measured
  # in WHOLE BLOOD by LC-MS/MS (Methods, "Tacrolimus dosing"), not plasma;
  # tacrolimus partitions extensively into erythrocytes, which is why
  # hematocrit enters the CL/F model at all.
  compartmentData <- list(
    depot   = list(analyte = "tacrolimus", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "tacrolimus", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list(
    POD = list(
      description        = "Postoperative day -- days elapsed since lung transplantation",
      units              = "days",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying within subject; the strongest single predictor of tacrolimus CL/F in Miano 2024",
        "(univariate Delta-OFV = -979.01, Table 3). Enters CL/F as an UN-normalised power term,",
        "POD^0.67, exactly as printed in the Results display equation -- unlike the hematocrit and",
        "weight terms in the same equation, POD carries no median divisor. The un-normalised reading",
        "is corroborated internally: at the dataset's median postoperative day (~7) the final model",
        "gives 3.69 * 7^0.67 = 13.6 L/h, which reproduces the covariate-free base-model CL/F of",
        "14.06 L/h (Table 2); a median-normalised reading would give 3.69 L/h, a 4-fold mismatch.",
        "The reference level of CL/F therefore corresponds to POD = 1 day. NOTE: the term is",
        "singular at POD = 0 (CL/F would be 0). This is not a practical restriction for the source",
        "dataset -- tacrolimus was started 12 to 24 hours postoperatively and troughs were drawn",
        "before the morning dose -- but datasets passed to this model must use POD >= 1.",
        "Miano 2024 collected data over postoperative days 0-14; median follow-up 13 days (IQR 11-14)."
      ),
      source_name        = "postoperative day"
    ),
    HCT = list(
      description        = "Hematocrit -- packed red blood cell volume fraction",
      units              = "%",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying within subject (Miano 2024 Methods: 'hematocrit (a marker for tacrolimus",
        "binding to red blood cells), modeled as a time-varying covariate'). Normalised power term",
        "(HCT/33)^-1.45 with the divisor 33 % equal to the derivation-cohort median hematocrit",
        "reported in Table 1 (33, IQR 30-38). Expressed as a PERCENT (0-100), matching the register",
        "canonical -- not as a fraction. The exponent is negative because tacrolimus partitions into",
        "erythrocytes: a higher red-cell volume sequesters more drug out of the plasma pool that is",
        "available for hepatic extraction, lowering apparent whole-blood clearance."
      ),
      source_name        = "hematocrit"
    ),
    TX_LUNG_BILAT = list(
      description        = "Bilateral (vs single) lung transplant indicator: 1 = bilateral / double lung graft, 0 = single / unilateral lung graft",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (single / unilateral lung transplant)",
      notes              = paste(
        "Time-fixed per subject. Every subject in the Miano 2024 cohort is a lung transplant",
        "recipient, so this covariate encodes graft LATERALITY, not organ identity (contrast the",
        "register's TX_LUNG, whose reference category is a non-lung organ). 66% of the derivation",
        "cohort and 68% of the validation cohort received bilateral grafts (Table 1). The effect on",
        "CL/F is TIME-VARYING: 0.48 on postoperative days 1-3 and 0.82 on days 4-14. The time-varying",
        "specification was required by the data -- a single time-invariant coefficient gave",
        "Delta-OFV = -2.87 (not significant) whereas the early/late split gave Delta-OFV = -576.71",
        "(Table 3). Miano 2024 interprets bilateral transplantation as a marker of a larger surgical",
        "insult (longer cardiopulmonary bypass, median sternotomy, greater ischemia-reperfusion",
        "risk) whose effect on clearance wanes as patients recover."
      ),
      source_name        = "transplant type"
    ),
    CYP3A5_EXPR = list(
      description        = "CYP3A5 expresser indicator: 1 if the patient carries at least one functional CYP3A5*1 allele (intermediate or extensive metabolizer), 0 if a poor metabolizer",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CYP3A5 poor metabolizer)",
      notes              = paste(
        "Time-fixed germline genotype. Miano 2024 genotyped three loss-of-function alleles --",
        "rs776746 (CYP3A5*3), rs10264272 (CYP3A5*6) and rs41303343 (CYP3A5*7) -- and assigned CPIC",
        "phenotypes: extensive metabolizer (*1/*1), intermediate metabolizer (*1/*3, *1/*6, *1/*7),",
        "poor metabolizer (*3/*3, *6/*6, *7/*7, *3/*6, *3/*7, *6/*7). Extensive and intermediate",
        "metabolizers were POOLED into a single group 'given small numbers in each category'",
        "(Methods, 'Genotypes'), which is exactly the expresser-equals-1 orientation of the register",
        "canonical. Note that this cohort's *1-carrier definition rests on three loss-of-function",
        "loci rather than rs776746 alone, so a dataset genotyped only at rs776746 will misclassify",
        "the small number of *6 and *7 carriers. 20% of the derivation cohort and 19% of the",
        "validation cohort were expressers (Table 1)."
      ),
      source_name        = "CYP3A5 metabolizer genotype"
    ),
    CONMED_VORICONAZOLE = list(
      description        = "Concomitant voriconazole exposure indicator (strong CYP3A inhibitor): 1 = treated, 0 = untreated",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no voriconazole)",
      notes              = paste(
        "Time-varying within subject and LAGGED BY 24 HOURS: Miano 2024 shifted the administration",
        "dates forward by one day 'to account for (1) the onset and offset of inhibition and (2) the",
        "timing of drug concentration monitoring', because troughs were drawn between 04:00 and",
        "06:00 and an inhibitor started later the same day would first affect the following",
        "morning's trough (Methods, 'Clinical covariates'). Datasets used with this model must apply",
        "the same 24-hour lag. Voriconazole is used for invasive-aspergillosis prophylaxis and was",
        "the single strongest co-medication effect (univariate Delta-OFV = -33.32; removal from the",
        "full model raised OFV by 338.01, Table 3). 35% of the derivation cohort were ever exposed."
      ),
      source_name        = "voriconazole"
    ),
    CONMED_AMIO = list(
      description        = "Concomitant amiodarone exposure indicator (CYP3A inhibitor): 1 = treated, 0 = untreated",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no amiodarone)",
      notes              = paste(
        "Time-varying within subject and lagged by 24 hours, on the same basis as",
        "CONMED_VORICONAZOLE (Miano 2024 Methods, 'Clinical covariates'). Amiodarone is used for",
        "prevention and treatment of postoperative atrial fibrillation in this population; 36% of",
        "the derivation cohort were ever exposed (Table 1)."
      ),
      source_name        = "amiodarone"
    ),
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Pretransplant body weight, time-fixed per subject. Normalised power terms (WT/75)^0.40 on",
        "CL/F and (WT/75)^0.79 on Vd/F; the divisor 75 kg equals the derivation-cohort median weight",
        "in Table 1 (75, IQR 59-88). Both exponents were ESTIMATED, not fixed at the allometric",
        "0.75 / 1.0 -- the CL/F exponent 0.40 (95% CI 0.11-0.69) and the Vd/F exponent 0.79",
        "(95% CI 0.41-1.18). NOTE a paper-internal inconsistency: the Results text quotes median",
        "weights of '77 kg vs 74 kg' for the two cohorts, whereas Table 1 reports 75 kg (derivation)",
        "and 77 kg (validation). The model uses 75 kg, the derivation-cohort median from Table 1,",
        "which is also the divisor printed in the Results display equation."
      ),
      source_name        = "weight"
    )
  )

  # Covariates that Miano 2024 SCREENED but did not retain in the final model.
  # Documented here so the paper's covariate search is auditable; these names
  # are deliberately absent from model(). Age and race were excluded before
  # screening because they were strongly collinear with cystic fibrosis and
  # CYP3A5 genotype respectively (Methods, "Clinical covariates").
  covariatesDataExcluded <- list(
    PGD = list(
      description        = "Primary graft dysfunction diagnosed by postoperative day 3 or earlier",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no primary graft dysfunction)",
      notes              = "Screened on CL/F both time-invariantly (Delta-OFV = -0.08) and with an early/late time-varying split (Delta-OFV = -141.26), and carried into the forward-selection sequence as Model 10, but dropped: adding it to Model 9 gave Delta-OFV = -3.50, below the 3.84 retention threshold, and Delta-AIC = +0.10 (Miano 2024 Table 3). 23% of the derivation cohort. Not retained in the final model; no point estimate is reported for it.",
      source_name        = "PGD"
    ),
    DIS_CF = list(
      description        = "Cystic fibrosis as the indication for transplantation",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-cystic-fibrosis indication)",
      notes              = "Screened univariately on CL/F with Delta-OFV = -0.10 (Miano 2024 Table 3) and not carried forward. Selected over age for screening because the two were strongly correlated and cystic fibrosis has the clearer biologic relationship with tacrolimus metabolism (Methods, 'Clinical covariates'). 11% of the derivation cohort.",
      source_name        = "cystic fibrosis"
    ),
    CONMED_FLUCONAZOLE = list(
      description        = "Concomitant fluconazole exposure indicator (CYP3A inhibitor)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no fluconazole)",
      notes              = "Screened univariately on CL/F with Delta-OFV = -10.02, and tested in forward selection as Model 8 (Delta-OFV = -3.06 versus Model 6, below the 3.84 retention threshold), so not retained (Miano 2024 Table 3). Time-varying and 24-hour-lagged in the source dataset, as for the other CYP inhibitors. 7% of the derivation cohort.",
      source_name        = "fluconazole"
    ),
    SNP_CYP3A4_RS35599367 = list(
      description        = "CYP3A4*22 (rs35599367, g.15389C>T in intron 6) reduced-function allele carrier indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no CYP3A4*22 allele)",
      notes              = "Screened univariately on CL/F with Delta-OFV = -7.54, and tested in forward selection as Model 7 (Delta-OFV = -1.36 versus Model 6), so not retained (Miano 2024 Table 3). The Discussion attributes the null multivariable result to the low number of carriers. Table 1 reports 29 (7%) of 270 derivation-cohort subjects as carriers -- note that 29/270 is 10.7%, not 7%, an unresolved arithmetic inconsistency in the source table.",
      source_name        = "CYP3A4 metabolizer genotype"
    ),
    AGE = list(
      description        = "Age at transplantation",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Listed among the pretransplant health-status covariates considered, but excluded from screening because age and cystic fibrosis were strongly correlated and cystic fibrosis was selected instead (Miano 2024 Methods, 'Clinical covariates'). Derivation-cohort median 61 years (IQR 51-66). No point estimate is reported.",
      source_name        = "age"
    ),
    RACE_BLACK = list(
      description        = "African American race indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (White or Other)",
      notes              = "Listed among the pretransplant health-status covariates considered, but excluded from screening because race and CYP3A5 genotype were strongly correlated and CYP3A5 genotype was selected instead (Miano 2024 Methods, 'Clinical covariates'). Derivation cohort: White 237 (88%), African American 23 (9%), Other 10 (3%). No point estimate is reported.",
      source_name        = "race"
    ),
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Listed among the pretransplant health-status covariates considered in Miano 2024 Methods, 'Clinical covariates', but no univariate screening result is reported for it in Table 3 and it does not appear in the final model. 40% of the derivation cohort were female. No point estimate is reported.",
      source_name        = "sex"
    )
  )

  population <- list(
    species           = "human",
    n_subjects        = 270L,
    n_studies         = 1L,
    n_observations    = 3143L,
    age_range         = "51-66 years (IQR); minimum 18 years by inclusion criterion",
    age_median        = "61 years",
    weight_range      = "59-88 kg (IQR)",
    weight_median     = "75 kg",
    sex_female_pct    = 39.6,
    race_ethnicity    = c(White = 88, `African American` = 9, Other = 3),
    disease_state     = "Adult lung transplant recipients during the first 14 postoperative days, on protocol immunosuppression. Approximately two-thirds received bilateral grafts; 23% developed primary graft dysfunction by postoperative day 3; 11% were transplanted for cystic fibrosis. All patients received corticosteroids (median 37 mg/day, IQR 32-44).",
    dose_range        = "Oral or sublingual tacrolimus twice daily (06:00 and 18:00), started 12-24 hours postoperatively, with doses at the discretion of the treating clinician and titrated to a whole-blood trough target of 8-12 ng/mL. Median initial dose 2 mg (IQR 2-2), Table 1. Patients treated with intravenous tacrolimus or with cyclosporine were excluded.",
    regions           = "United States (single centre; Hospital of the University of Pennsylvania, Philadelphia).",
    cyp3a5_distribution = "CYP3A5 intermediate or extensive metabolizer (at least one functional *1 allele) 53/270 (20%) in the derivation cohort and 22/114 (19%) in the validation cohort; the remainder poor metabolizers. Phenotypes assigned per CPIC from rs776746 (*3), rs10264272 (*6) and rs41303343 (*7).",
    hematocrit_median = "33% (IQR 30-38); time-varying",
    validation_cohort = "An additional 114 patients (1,279 concentrations) were held out for external validation and were NOT used to estimate the parameters in this file. Mean prediction error in the validation cohort was 36.4% (95% CI 30.8-41.9) and median prediction error 7.2% (IQR -29.3 to 70.53); 34.7% and 59.1% of population-predicted concentrations fell within 2 and 4 ng/mL of the observed value respectively.",
    notes             = paste(
      "Retrospective population PK analysis of patients enrolled in the multicentre Lung Transplant",
      "Outcomes Group (LTOG) cohort at the University of Pennsylvania between November 2008 and",
      "August 2018, merged with electronic health record data. 384 genotyped patients with at least",
      "one tacrolimus concentration were randomly split approximately 70/30 into a derivation cohort",
      "(n = 270, the parameters in this file) and a validation cohort (n = 114). Data were collected",
      "over the first 14 postoperative days; median follow-up 13 days (IQR 11-14) with a median of",
      "13 samples per patient (IQR 10-14). Sampling was almost entirely TROUGH-ONLY -- median time",
      "after dose 11.1 hours (IQR 10.1-11.9) -- which is why the authors fixed ka a priori and",
      "specified a one-compartment model rather than estimating absorption or distribution.",
      "Tacrolimus was measured in whole blood by LC-MS/MS with an analytical range of 1-30 ng/mL",
      "(within- and between-run CV < 9.2%). NONMEM 7.5.1, FOCE with interaction. Baseline",
      "demographics are in Table 1; parameter estimates in Table 2; the covariate search in Table 3."
    )
  )

  ini({
    # Structural PK -- Miano 2024 Table 2, "Final model" column. Time in hours;
    # apparent clearance CL/F in L/h; apparent volume Vd/F in L. Both are
    # reported relative to bioavailability (Methods, "Base model development":
    # "Clearance (CL) and volume of distribution (Vd) are reported relative to
    # bioavailability (F)"), so F is absorbed into the parameters and no
    # f(depot) is assigned in model().
    #
    # The reference subject for lcl is a single-lung recipient on postoperative
    # day 1 with hematocrit 33%, CYP3A5 poor-metabolizer genotype, no
    # voriconazole, no amiodarone, and weight 75 kg -- i.e. every covariate
    # factor in the Results display equation equal to 1.
    lka <- fixed(log(4.5)) ; label("Absorption rate constant ka (1/h)")                                        # Miano 2024 Table 2: "Ka, liter/hour | 4.5 Fixed" -- fixed a priori from reference 16, NOT estimated. The table's unit string "liter/hour" is a typographical error: ka is a first-order rate constant and must have units of 1/h. A sensitivity analysis with ka fixed to 0.58 1/h (Table S1) gave "negligible effects on model parameters" (Discussion).
    lcl <- log(3.69)       ; label("Apparent oral clearance CL/F at POD 1 d, HCT 33%, single-lung graft, CYP3A5 poor metabolizer, no voriconazole, no amiodarone, WT 75 kg (L/h)") # Miano 2024 Table 2 final model CL/F = 3.69 L/h (95% CI 2.93-4.45; bootstrap median 3.68, 2.5-97.5 pct 2.77-4.47)
    lvc <- log(642)        ; label("Apparent volume of distribution Vd/F at WT 75 kg (L)")                     # Miano 2024 Table 2 final model Vd/F = 642 L (95% CI 575-709; bootstrap median 644, 2.5-97.5 pct 604-779)

    # Covariate effects on CL/F -- Miano 2024 Results display equation:
    #   CL/F (L/h) = 3.69 * POD^0.67 * (HCT/33)^-1.45
    #                * [1.0 single | 0.48 bilateral days 1-3 | 0.82 bilateral days 4-14]
    #                * [1.0 CYP3A5 poor metabolizer | 1.78 intermediate or extensive]
    #                * [1.0 voriconazole untreated | 0.39 voriconazole treated]
    #                * [1.0 amiodarone untreated | 0.77 amiodarone treated]
    #                * (WT/75)^0.40
    # The binary factors follow the multiplicative form the paper specifies in
    # Methods, "Modeling procedures": TVP = theta_TVP * (theta_cov)^COV with
    # covariates coded 1 = present, 0 = absent. POD, HCT and WT follow the
    # normalised power form TVP = theta_TVP * (cov/cov_median)^theta_cov from
    # the same paragraph -- except that POD is printed WITHOUT a median
    # divisor; see covariateData$POD$notes for the internal check confirming
    # that the un-normalised form is the intended one.
    e_pod_cl                  <- 0.67  ; label("Power exponent on postoperative day for CL/F (unitless)")                                # Miano 2024 Table 2 "Postoperative day" = 0.67 (95% CI 0.58-0.76)
    e_hct_cl                  <- -1.45 ; label("Power exponent on (HCT / 33 %) for CL/F (unitless)")                                     # Miano 2024 Table 2 "Hematocrit, %" = -1.45 (95% CI -1.84 to -1.06)
    e_tx_lung_bilat_cl_pod13  <- 0.48  ; label("Bilateral (vs single) lung graft multiplicative factor on CL/F, postoperative days 1-3")  # Miano 2024 Table 2 "Transplant type / Days 1-3" = 0.48 (95% CI 0.35-0.61)
    e_tx_lung_bilat_cl_pod414 <- 0.82  ; label("Bilateral (vs single) lung graft multiplicative factor on CL/F, postoperative days 4-14") # Miano 2024 Table 2 "Transplant type / Days 4-14" = 0.82 (95% CI 0.70-0.94)
    e_cyp3a5_expr_cl          <- 1.78  ; label("CYP3A5 intermediate/extensive metabolizer multiplicative factor on CL/F")                # Miano 2024 Table 2 "CYP3A5 genotype" = 1.78 (95% CI 1.44-2.12)
    e_conmed_voriconazole_cl  <- 0.39  ; label("Concomitant voriconazole multiplicative factor on CL/F")                                 # Miano 2024 Table 2 "Voriconazole exposure" = 0.39 (95% CI 0.32-0.47)
    e_conmed_amio_cl          <- 0.77  ; label("Concomitant amiodarone multiplicative factor on CL/F")                                   # Miano 2024 Table 2 "Amiodarone exposure" = 0.77 (95% CI 0.69-0.85)
    e_wt_cl                   <- 0.40  ; label("Power exponent on (WT / 75 kg) for CL/F (unitless)")                                     # Miano 2024 Table 2 "Covariates on CL / Weight, kg" = 0.40 (95% CI 0.11-0.69) -- estimated, not fixed at 0.75

    # Covariate effect on Vd/F -- Miano 2024 Results display equation:
    #   Vd (L) = 642 * (WT/75)^0.79
    e_wt_vc                   <- 0.79  ; label("Power exponent on (WT / 75 kg) for Vd/F (unitless)")                                     # Miano 2024 Table 2 "Covariates on Vd / Weight, kg" = 0.79 (95% CI 0.41-1.18) -- estimated, not fixed at 1.0

    # Inter-individual variability -- Miano 2024 Table 2 "Interindividual
    # variability" rows, final-model column, taken as NONMEM OMEGA VARIANCES on
    # the log scale for the exponential model the paper specifies (Methods,
    # "Base model development": "We described random interindividual
    # variability using an exponential variance model"). The table prints these
    # two rows as bare decimals with NO scale label, unlike the residual rows
    # immediately below them which are explicitly labelled "%" and "ng/ml"; the
    # unlabelled reading is the raw NONMEM element, which is a variance.
    #   CL/F  omega^2 = 0.29 -> CV = sqrt(exp(0.29) - 1) = 58%
    #   Vd/F  omega^2 = 0.54 -> CV = sqrt(exp(0.54) - 1) = 85%
    # Operator-ratified reading (sidecar request 001, question 1). See the
    # vignette's "Assumptions and deviations" section, which records the
    # alternative SD reading (CV 30% and 58%) so a reader holding the control
    # stream can settle it.
    #
    # Methods state the random effects were "assumed to follow a multivariate
    # normal covariance structure", i.e. a full OMEGA BLOCK, but Table 2 reports
    # only the two diagonal elements and no covariance or correlation anywhere
    # in the paper. The off-diagonal is therefore encoded as fixed(0) rather
    # than invented -- this is a KNOWN DEVIATION from the published structure,
    # not a claim that the authors fixed it to zero.
    etalcl + etalvc ~ c(0.29,
                        fixed(0), 0.54)                                                                                                 # Miano 2024 Table 2 IIV CL/F = 0.29 (95% CI 0.23-0.34, shrinkage 4.8%), IIV Vd/F = 0.54 (95% CI 0.43-0.65, shrinkage 7.6%); off-diagonal not reported

    # Residual unexplained variability -- Miano 2024 Table 2 "Residual
    # variability" rows, final-model column. Both are on the SD scale: the
    # proportional row is labelled "%" and the additive row "ng/ml" (not
    # "(ng/ml)^2"). The base-model row proves the convention -- its point
    # estimate is printed as 13.01 while its own confidence interval is printed
    # as (0.09, 0.17), i.e. the underlying NONMEM number was 0.1301 and only the
    # point estimate was multiplied by 100. So 7.73 means a fractional
    # proportional SD of 0.0773. Corroboration: the LC-MS/MS assay's within- and
    # between-run CVs are reported as < 9.2% (Methods), and 7.73% sits just
    # under that. Operator-ratified reading (sidecar request 001, question 3).
    propSd <- 0.0773 ; label("Proportional residual error (fraction)")   # Miano 2024 Table 2 final model "Proportional, %" = 7.73 (95% CI 6.23-9.23, shrinkage 7.5%)
    addSd  <- 1.80   ; label("Additive residual error (ng/mL)")          # Miano 2024 Table 2 final model "Additive, ng/ml" = 1.80 (95% CI 1.08-2.52, shrinkage 14.4%)
  })

  model({
    # 1. Derived covariate terms.
    #
    # The bilateral-graft effect on CL/F is time-varying: Miano 2024 estimated
    # separate coefficients for postoperative days 1-3 and days 4-14 (Table 2
    # and the Results display equation). The Methods describe the same split as
    # "days 0-3 and 4-14"; either wording places the breakpoint between day 3
    # and day 4, which is where the authors' prior work found an inflection in
    # the tacrolimus concentration:dose ratio.
    txLungBilatEffect <- e_tx_lung_bilat_cl_pod414
    if (POD < 4) {
      txLungBilatEffect <- e_tx_lung_bilat_cl_pod13
    }

    # 2. Individual PK parameters.
    #
    # POD enters as an un-normalised power term, reproducing the Results display
    # equation verbatim. It is singular at POD = 0; see covariateData$POD$notes.
    cl <- exp(lcl + etalcl) *
      POD^e_pod_cl *
      (HCT / 33)^e_hct_cl *
      txLungBilatEffect^TX_LUNG_BILAT *
      e_cyp3a5_expr_cl^CYP3A5_EXPR *
      e_conmed_voriconazole_cl^CONMED_VORICONAZOLE *
      e_conmed_amio_cl^CONMED_AMIO *
      (WT / 75)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 75)^e_wt_vc
    ka <- exp(lka)

    # 3. Micro-constants.
    kel <- cl / vc

    # 4. One-compartment oral / sublingual PK with first-order absorption. The
    # authors specified one compartment a priori because the dataset is almost
    # entirely trough concentrations, which "precluded estimation of
    # intercompartmental transfer" (Methods, "Base model development").
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # 5. Bioavailability is not assigned: CL/F and Vd/F already carry F.

    # 6. Observation and error. Doses are in mg and vc is in L, so central/vc is
    # in mg/L; multiply by 1000 to reach the ng/mL of the whole-blood assay.
    Cc <- central / vc * 1000
    Cc ~ add(addSd) + prop(propSd)
  })
}
