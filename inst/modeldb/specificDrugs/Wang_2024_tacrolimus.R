Wang_2024_tacrolimus <- function() {
  description <- paste0(
    "One-compartment population PK model with first-order absorption for ",
    "twice-daily oral tacrolimus (Prograf) trough concentrations in Chinese ",
    "adult renal transplant recipients (Wang 2024). The absorption rate ",
    "constant ka is held at 3.86 1/h because only pre-dose trough ",
    "concentrations (C0) were available. Apparent oral clearance CL/F carries ",
    "two exponential covariate effects: the CYP3A5 *3-allele count (0 for ",
    "*1/*1, 1 for *1/*3, 2 for *3/*3; coefficient -0.348, so CL/F falls to ",
    "70.6% and 49.9% of the *1/*1 value for *1/*3 and *3/*3 respectively) and ",
    "hematocrit expressed as a volume fraction (coefficient -0.122). Apparent ",
    "central volume V/F carries no covariate. Inter-individual variability is ",
    "exponential and diagonal on CL/F and V/F. Residual error is additive on ",
    "log-transformed concentrations, which is proportional error on the ",
    "linear concentration scale. Wang 2024 also builds MLP, SVR, and XGBoost ",
    "machine-learning predictors on top of this model's post hoc individual ",
    "predictions; those are not ODE models and are outside the scope of this ",
    "model file."
  )
  reference <- paste0(
    "Wang YP, Lu XL, Shao K, Shi HQ, Zhou PJ, Chen B. Improving prediction of ",
    "tacrolimus concentration using a combination of population ",
    "pharmacokinetic modeling and machine learning in chinese renal ",
    "transplant recipients. Front Pharmacol. 2024;15:1389271. ",
    "doi:10.3389/fphar.2024.1389271."
  )
  vignette <- "Wang_2024_tacrolimus"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    HCT = list(
      description        = "Hematocrit -- packed red blood cell volume fraction",
      units              = "L/L (volume fraction)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Time-varying; recorded on each day of tacrolimus therapeutic drug ",
        "monitoring (Wang 2024 Methods 2.1). Wang 2024 Table 1 reports HCT as ",
        "a volume FRACTION (training set 0.29 +/- 0.056, test set 0.27 +/- ",
        "0.051), NOT as percent, and the printed final-model equation ",
        "CL/F = 70.6 * exp(CYP3A5 * -0.348) * exp(HCT * -0.122) consumes that ",
        "fraction directly. The canonical-register HCT entry's units (%) are ",
        "therefore explicitly overridden here, following the same override ",
        "already made in Andrews_2017_tacrolimus.R. To use a dataset that ",
        "records HCT in percent, multiply the column by 0.01 before passing ",
        "it to this model. The fraction-versus-percent reading is not a ",
        "judgement call: substituting the training-set mean covariates into ",
        "the final-model equation reproduces the covariate-free structural ",
        "model's CL/F of 41.1 L/h only under the fraction reading ",
        "(70.6 * exp(-0.122 * 0.29) * 0.618 = 42.1 L/h, 2.4% off), whereas ",
        "the percent reading gives 1.27 L/h (97% off). See the vignette's ",
        "Assumptions and deviations section. Higher hematocrit means more ",
        "erythrocyte binding of tacrolimus (tacrolimus is extensively ",
        "red-cell bound) and hence lower apparent whole-blood clearance, ",
        "which is why the coefficient is negative. Wang 2024 Discussion ",
        "reports that adding HCT decreased the inter-individual variability ",
        "in CL/F by 6.21%."
      ),
      source_name        = "HCT"
    ),
    CYP3A5_STAR1_HET = list(
      description        = "CYP3A5*1/*3 heterozygote indicator (one functional CYP3A5*1 allele)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not a *1/*3 heterozygote)",
      notes              = paste0(
        "Time-fixed per subject (germline genotype at rs776746, determined by ",
        "PCR-sequencing; Wang 2024 Methods 2.2). 1 = subject is CYP3A5*1/*3; ",
        "0 = otherwise (the union of *1/*1 homozygotes and *3/*3 ",
        "nonexpressers). Paired with CYP3A5_STAR1_HOM. Wang 2024 does NOT ",
        "fit an independent coefficient per genotype stratum: it assigns the ",
        "ordinal scores 0, 1, 2 to *1/*1, *1/*3, *3/*3 (Wang 2024 Methods 2.3 ",
        "and Results 3.2) and enters that score linearly into ",
        "exp(CYP3A5 * -0.348). The score is the count of nonfunctional *3 ",
        "alleles, which model() reconstructs exactly from this pair of ",
        "registered binary indicators as ",
        "2 - (2 * CYP3A5_STAR1_HOM + CYP3A5_STAR1_HET). Encoding via the ",
        "existing paired binaries (the Passey_2011_tacrolimus.R precedent) ",
        "avoids registering a redundant CYP3A5 genotype-count canonical, ",
        "which inst/references/covariate-columns.md explicitly discourages ",
        "under the CYP3A5_EXPR entry. Wang 2024 Table 1: 35/103 (34%) of the ",
        "training set and 8/24 (33.3%) of the test set were *1/*3."
      ),
      source_name        = "CYP3A5 *1/*3"
    ),
    CYP3A5_STAR1_HOM = list(
      description        = "CYP3A5*1/*1 homozygote indicator (two functional CYP3A5*1 alleles)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not a *1/*1 homozygote)",
      notes              = paste0(
        "Time-fixed per subject (germline genotype at rs776746). 1 = subject ",
        "is CYP3A5*1/*1; 0 = otherwise (the union of *1/*3 heterozygotes and ",
        "*3/*3 nonexpressers). Paired with CYP3A5_STAR1_HET; see that entry ",
        "for how the two indicators reconstruct Wang 2024's ordinal 0/1/2 ",
        "genotype score. CYP3A5_STAR1_HET = 0 AND CYP3A5_STAR1_HOM = 0 is the ",
        "*3/*3 nonexpresser group, which carries the score 2 -- note that ",
        "this is the OPPOSITE orientation from the Passey_2011_tacrolimus.R ",
        "use of the same pair, where *3/*3 is the reference with a factor of ",
        "1.0. Wang 2024 anchors its equation at *1/*1 (score 0) instead, so ",
        "the reported CL/F of 70.6 L/h is the *1/*1 value. Wang 2024 Table 1: ",
        "10/103 (9.7%) of the training set and 2/24 (8.3%) of the test set ",
        "were *1/*1."
      ),
      source_name        = "CYP3A5 *1/*1"
    )
  )

  # Covariates that Wang 2024 collected and screened during forward
  # inclusion / backward elimination but did NOT retain in the final PPK
  # model. Several of them (notably POD) rank highly in the paper's SHAP
  # analysis of the machine-learning layer, but the paper is explicit that
  # "only CYP3A5 genotype and HCT were proved to be the covariates of TAC PPK
  # model" (Wang 2024 Discussion). They are documented here so the provenance
  # of the covariate screen is preserved without declaring covariates that
  # model() never references. Platelet count and total bilirubin were also
  # screened but have no canonical entry in
  # inst/references/covariate-columns.md, so they are recorded in
  # population$notes instead of here.
  covariatesDataExcluded <- list(
    POD = list(
      description = "Postoperative date -- days between the transplant operation and the day of data collection",
      units       = "days",
      type        = "continuous",
      notes       = "Screened as a covariate on CL/F and V/F but not retained in the final model. Wang 2024 Table 1: 130 +/- 231 days (training), 175 +/- 242 days (test); overall range 3-1622 days. Ranked second (behind the PPK basic model's post hoc IPRE) in the SHAP importance analysis of the machine-learning layer."
    ),
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened but not retained. Wang 2024 Table 1: 63.3 +/- 12.9 kg (training), 62.6 +/- 11.0 kg (test). No allometric scaling appears in the final model."
    ),
    AGE = list(
      description = "Recipient age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened but not retained. Wang 2024 Table 1: 41.2 +/- 11.2 years (training), 44.7 +/- 8.91 years (test)."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened but not retained. Wang 2024 Methods 2.3 codes sex as 0 for male and 1 for female, which matches the canonical SEXF orientation. Table 1: 39/103 female (training), 7/24 female (test)."
    ),
    RBC = list(
      description = "Red blood cell count",
      units       = "10^12 cells/L",
      type        = "continuous",
      notes       = "Screened but not retained. Wang 2024 Table 1: 3.06 +/- 0.61 (training), 2.98 +/- 0.51 (test). Correlated with the retained HCT covariate."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened but not retained. Wang 2024 Table 1: 34.0 +/- 3.69 g/L (training), 33.1 +/- 4.04 g/L (test)."
    ),
    ALP = list(
      description = "Alkaline phosphatase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened but not retained. Wang 2024 Table 1 reports ALP with the units 'mmol/L', which is not a valid unit for an enzyme activity; U/L is the standard reporting unit and the tabulated magnitudes (58.2 +/- 24.9 training, 63.0 +/- 44.6 test) are consistent with U/L."
    ),
    CRCL = list(
      description = "Creatinine clearance",
      units       = "mL/min/1.73 m^2",
      type        = "continuous",
      notes       = "Screened but not retained. Wang 2024 Methods 2.1 lists creatinine clearance among the collected clinical data; Table 1 tabulates serum creatinine (212.5 +/- 208.5 umol/L training, 305.5 +/- 340.6 umol/L test) rather than the derived clearance."
    )
  )

  compartmentData <- list(
    depot   = list(analyte = "tacrolimus", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "tacrolimus", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  population <- list(
    species          = "human",
    n_subjects       = 103L,
    n_studies        = 1L,
    n_concentrations = 2041L,
    age_range        = "adults; Table 1 mean +/- SD 41.2 +/- 11.2 years (training set), 44.7 +/- 8.91 years (test set); 42.2 +/- 11.0 years across all 127 recipients",
    weight_range     = "Table 1 mean +/- SD 63.3 +/- 12.9 kg (training set), 62.6 +/- 11.0 kg (test set); 62.6 +/- 12.2 kg across all 127 recipients",
    sex_female_pct   = 37.9,
    race_ethnicity   = c(Asian = 100),
    disease_state    = paste0(
      "Adult Chinese recipients of a first renal transplant meeting standard ",
      "renal donor criteria, on a triple immunosuppressive regimen of ",
      "tacrolimus + mycophenolate mofetil + corticosteroids. Excluded: ",
      "combined organ transplantation, panel-reactive-antibody positivity, ",
      "tacrolimus allergy or intolerance, pregnancy or lactation. Concomitant ",
      "medications in the training set: calcium antagonists 81 (78.6%), ",
      "proton pump inhibitors 102 (99%), voriconazole 21 (20.4%)."
    ),
    dose_range       = paste0(
      "Oral tacrolimus (Prograf, Astellas) started at 0.1 mg/kg/day given as ",
      "two divided doses q12h, then titrated to a trough target of 10-13 ",
      "ng/mL during the first month post-transplant and 5-9 ng/mL ",
      "thereafter. Tacrolimus was administered between 3 and 1622 days after ",
      "transplantation."
    ),
    regions          = "China (single centre: Ruijin Hospital, Shanghai Jiaotong University School of Medicine)",
    cyp3a5_genotype  = c(`*1/*1` = 9.7, `*1/*3` = 34.0, `*3/*3` = 56.3),
    sampling_design  = paste0(
      "Routine therapeutic drug monitoring only: a single pre-dose trough ",
      "sample (C0) drawn at 08:00 immediately before the morning dose. No ",
      "post-dose or rich sampling was collected, which is why a ",
      "one-compartment structural model was used and ka was held constant ",
      "rather than estimated (Wang 2024 Results 3.2). Whole-blood tacrolimus ",
      "was assayed by enzyme-multiplied immunoassay (Syva VivaEmit 2000, ",
      "Siemens), quantification range 2-50 ng/mL."
    ),
    notes            = paste0(
      "2041 trough concentrations from 127 recipients were split at random ",
      "into a training set (103 patients, 80%) and a test set (24 patients, ",
      "20%). The parameter estimates encoded here are the 'Final model' ",
      "column of Wang 2024 Table 2, which was fitted to the TRAINING set; ",
      "population$n_subjects is therefore 103, not 127. Table 2 also reports ",
      "a 'Test set final model' column (V/F 2330 L, CL/F 114 L/h, HCT ",
      "coefficient -0.161, CYP3A5 coefficient -0.395) obtained by re-fitting ",
      "the same structure to the 24-patient test set; that refit is a ",
      "validation exercise, not the paper's final model, and is not encoded ",
      "here. Software: NONMEM 7.4.1 (FOCE), Pirana 23.1.1, PsN 5.3.1, R ",
      "4.3.0. Covariate selection used forward inclusion (dOFV <= -6.63, ",
      "p < 0.01, df = 1) and backward elimination (dOFV >= 7.88, p < 0.005, ",
      "df = 1). Additional covariates that were collected and screened but ",
      "have no canonical column name in inst/references/covariate-columns.md ",
      "are platelet count (150.6 +/- 54.3 x 10^9/L training) and total ",
      "bilirubin (11.3 +/- 4.49 umol/L training); neither was retained. ",
      "White blood cell count and blood urea nitrogen were also tabulated. ",
      "Note that Table 1's test-set BUN entry (243.0 +/- 14.4 mmol/L) is ",
      "implausible for that unit and is presumably a typographical error in ",
      "the source; BUN is not used by this model."
    )
  )

  ini({
    # ----- Structural parameters (Wang 2024 Table 2, "Final model" column) -----

    # Absorption rate constant. Held constant during estimation because only
    # pre-dose trough concentrations were available, so ka is not
    # identifiable from the data. NOTE a source-internal inconsistency: Wang
    # 2024 Table 2 reports 3.86 1/h in all three model columns (theta1), while
    # Results 3.2 states "The value of ka was fixed at 3.84 h-1". Table 2 is
    # the paper's designated parameter table and states 3.86 three times, so
    # 3.86 is used here; the 0.5% difference has no material effect on a
    # trough-only model. See the vignette's Assumptions and deviations
    # section.
    lka <- fixed(log(3.86)); label("Absorption rate constant ka (1/h)")  # Wang 2024 Table 2, theta1 = 3.86 (fixed)

    # Apparent central volume of distribution. No covariate was retained on
    # V/F. Structural (covariate-free) model reported 2620 L (RSE 10.3%).
    lvc <- log(2560); label("Apparent central volume of distribution V/F (L)")  # Wang 2024 Table 2, final model theta2 = 2560 L, RSE 10.7%

    # Apparent oral clearance for the reference subject. Wang 2024's
    # final-model equation is uncentred in both covariates, so this value is
    # CL/F at CYP3A5 score 0 (the *1/*1 genotype) AND hematocrit 0 L/L. The
    # hematocrit anchor is non-physiological, but it is what the printed
    # equation says; at the training-set mean HCT of 0.29 L/L the exponential
    # factor is exp(-0.122 * 0.29) = 0.965, so the *1/*1 typical CL/F at mean
    # hematocrit is 68.1 L/h. The structural model without covariates reported
    # CL/F = 41.1 L/h (RSE 14.2%), which the final model reproduces at the
    # training-set mean covariates (see the HCT covariateData notes).
    lcl <- log(70.6); label("Apparent oral clearance CL/F at CYP3A5*1/*1 and hematocrit 0 L/L (L/h)")  # Wang 2024 Table 2, final model theta3 = 70.6 L/h, RSE 9.75%

    # ----- Covariate effects on CL/F (Wang 2024 Results 3.2 equation) -----
    # Printed equation:
    #   CL/F = 70.6 * e^(CYP3A5 * -0.348) * e^(HCT * -0.122)
    # where CYP3A5 is the ordinal score 0 / 1 / 2 for *1/*1 / *1/*3 / *3/*3
    # (Wang 2024 Methods 2.3 and Results 3.2) and HCT is the hematocrit
    # volume fraction (Table 1 reports 0.29, not 29).

    # The exponential coefficient reproduces the genotype effect sizes stated
    # in the Discussion exactly: exp(-0.348) = 0.706 and exp(-0.348 * 2) =
    # 0.4985, matching "the CL/F of CYP3A5*1/*3 and *3/*3 patients were 70.6%
    # and 49.9% of those with the *1/*1 genotype".
    e_cyp3a5_cl <- -0.348; label("Exponential coefficient of the CYP3A5*3-allele count (0/1/2) on CL/F (unitless)")  # Wang 2024 Table 2, final model theta5 = -0.348, RSE 12.2%

    e_hct_cl <- -0.122; label("Exponential coefficient of hematocrit (L/L) on CL/F (unitless)")  # Wang 2024 Table 2, final model theta4 = -0.122, RSE 47.1%

    # ----- Inter-individual variability (Wang 2024 Table 2, "Final model") -----
    # IIV is exponential: P_i = TV(P_i) * exp(eta_i), eta_i ~ N(0, omega^2)
    # (Wang 2024 Methods 2.3, equation 1). Table 2 labels the two IIV rows
    # "omega V/F (%)" and "omega CL/F (%)", i.e. it reports omega itself --
    # the standard deviation on the log scale -- expressed as a percentage.
    # The variance entered below is therefore (percentage / 100)^2 and NOT
    # the log-normal CV% back-transform log(1 + CV^2). The Discussion
    # corroborates the SD reading: it states that adding CYP3A5 genotype
    # produced "an 8.85% decrease in interindividual variation in CL/F", and
    # 31.8 - 8.85 = 22.95, which is the reported final omega CL/F of 23.0 --
    # an arithmetic difference of percentage points on the omega scale. The
    # alternative CV% reading is discussed in the vignette's Assumptions and
    # deviations section.
    etalvc ~ 0.4225   # Wang 2024 Table 2, final model eta1: omega V/F = 65.0%, RSE 22.5%; variance = 0.650^2
    etalcl ~ 0.0529   # Wang 2024 Table 2, final model eta2: omega CL/F = 23.0%, RSE 23.1%; variance = 0.230^2

    # ----- Residual error (Wang 2024 Table 2, "Final model") -----
    # Wang 2024 Methods 2.3, equation 2: ln(Cobs) = ln(Cpred) + epsilon, i.e.
    # additive residual error on the log-transformed concentration scale.
    # Additive-on-log is proportional error in nlmixr2's linear concentration
    # space, so the 35.6% of Table 2 becomes propSd = 0.356.
    propSd <- 0.356; label("Proportional residual error (fraction)")  # Wang 2024 Table 2, final model delta = 35.6%, RSE 12.0%
  })

  model({
    # ----- 1. Derived covariate terms -----
    # Wang 2024 scores the CYP3A5 genotype as the count of nonfunctional *3
    # alleles: *1/*1 = 0, *1/*3 = 1, *3/*3 = 2. The registered binary
    # indicators flag the two functional-allele strata, so the number of
    # functional *1 alleles is (2 * homozygote + 1 * heterozygote) and the
    # paper's score is 2 minus that. Check: *1/*1 -> HOM 1, HET 0 -> 2-2 = 0;
    # *1/*3 -> HOM 0, HET 1 -> 2-1 = 1; *3/*3 -> HOM 0, HET 0 -> 2-0 = 2.
    cyp3a5Star1Count <- 2 * CYP3A5_STAR1_HOM + CYP3A5_STAR1_HET
    cyp3a5Star3Count <- 2 - cyp3a5Star1Count

    # ----- 2. Individual PK parameters -----
    ka <- exp(lka)

    # CL/F = 70.6 * exp(CYP3A5 * -0.348) * exp(HCT * -0.122) * exp(eta_CL).
    # HCT is consumed as the volume fraction (~0.29), not as percent -- see
    # the HCT covariateData notes for the arithmetic that settles this.
    cl <- exp(lcl + etalcl) *
      exp(e_cyp3a5_cl * cyp3a5Star3Count) *
      exp(e_hct_cl * HCT)

    vc <- exp(lvc + etalvc)

    # ----- 3. Micro-constants -----
    kel <- cl / vc

    # ----- 4. ODE system (one compartment, first-order absorption) -----
    # Bioavailability is not identifiable from oral-only data, so Wang 2024
    # reports the apparent parameters CL/F and V/F and no F term is encoded.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # ----- 5. Observation and error -----
    # Doses are in mg and V/F is in L, so central / vc is mg/L. Wang 2024
    # reports concentrations in ng/mL; 1 mg/L = 1 ug/mL = 1000 ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
