Lane_2011_warfarin_s <- function() {
  description <- "S-warfarin population PK (1-compartment, first-order absorption) in adults on long-term warfarin therapy (Lane 2011). Bodyweight, age, sex, and CYP2C9 diplotype influence apparent clearance; volume of distribution carries no covariates. Block correlation between random effects on CL and V. R-warfarin is reported separately in the same paper (modellib('Lane_2011_warfarin_r'))."
  reference <- paste(
    "Lane S, Al-Zubiedi S, Hatch E, Matthews I, Jorgensen AL, Deloukas P, Daly AK,",
    "Park BK, Aarons L, Ogungbenro K, Kamali F, Hughes D, Pirmohamed M.",
    "The population pharmacokinetics of R- and S-warfarin: effect of genetic and",
    "clinical factors.",
    "Br J Clin Pharmacol. 2012;73(1):66-76.",
    "doi:10.1111/j.1365-2125.2011.04051.x.",
    "PMID: 21692829.",
    "PK parameters and CYP2C9 / age / sex / weight effects from Table 3 (final",
    "covariate model); structural-equation form (CL_i = theta_CL * (WT/70)^theta_wgt",
    "* (1 + theta_age*(AGE-69.8)) * theta_CYP2C9 * theta_gender * exp(eta_CL)) from",
    "the Results paragraph following Tables 2 and 3."
  )
  vignette <- "Lane_2011_warfarin"

  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at baseline. Power scaling on CL with reference weight 70 kg (Lane 2011 final S-warfarin equation; the typical CL of 0.144 L/h is reported for a 70-kg, 69.8-year-old woman with CYP2C9 *1/*1). The Lane 2011 cohort had mean (range) body weight 80.7 (36-172) kg (Table 1).",
      source_name        = "WGT"
    ),
    AGE = list(
      description        = "Age at baseline",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at baseline. Linear fractional effect on CL with reference age 69.8 years (Lane 2011 final S-warfarin equation). The Lane 2011 cohort had mean (range) age 66.4 (19-95) years (Table 1); the 69.8 y reference value used in the equation is what the paper reports the typical CL at, not the cohort mean.",
      source_name        = "AGE"
    ),
    SEXF = list(
      description        = "Sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (female -- the reference category in the Lane 2011 final S-warfarin model; CL multiplier is 1.00 for women and 1.12 for men)",
      notes              = "Lane 2011 source coding was inverse (sex = 1 for man, 0 for woman; Methods 'Covariate selection and models'). Converted to canonical SEXF (1 = female, 0 = male). The CL multiplier for the non-reference (male) category is encoded as e_male_cl = 1.12 in ini() and applied as 1 + (e_male_cl - 1) * (1 - SEXF) in model(). Cohort sex distribution (Table 1): 58% male, 42% female.",
      source_name        = "SEX"
    ),
    CYP2C9_S1_COUNT = list(
      description        = "Count of CYP2C9*1 (wild-type) alleles per subject (0, 1, or 2)",
      units              = "(count, 0/1/2)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Paired with CYP2C9_S2_COUNT and CYP2C9_S3_COUNT; the three counts sum to 2 for any subject with known CYP2C9 genotype (CYP2C9_MISSING == 0). Cohort distribution (Table 1, n=306): *1/*1 63.7%, *1/*2 19.3%, *1/*3 9.5%, *2/*2 0.3%, *2/*3 2.0%, *3/*3 0.6%, missing 4.6%.",
      source_name        = "CYP2C9 (genotype string)"
    ),
    CYP2C9_S2_COUNT = list(
      description        = "Count of CYP2C9*2 reduced-function alleles per subject (0, 1, or 2)",
      units              = "(count, 0/1/2)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "See CYP2C9_S1_COUNT.",
      source_name        = "CYP2C9 (genotype string)"
    ),
    CYP2C9_S3_COUNT = list(
      description        = "Count of CYP2C9*3 reduced-function alleles per subject (0, 1, or 2)",
      units              = "(count, 0/1/2)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "See CYP2C9_S1_COUNT.",
      source_name        = "CYP2C9 (genotype string)"
    ),
    CYP2C9_MISSING = list(
      description        = "Binary indicator: CYP2C9 genotype not measured / not available",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CYP2C9 genotype known)",
      notes              = "Lane 2011 fits a separate CL multiplier (0.782) for the missing-genotype subgroup (14 of 306 S-warfarin subjects; Table 1 footnote). When CYP2C9_MISSING == 1 the three CYP2C9_S{1,2,3}_COUNT columns must all be 0 so the diplotype-indicator products evaluate to zero and the missing-multiplier is applied instead.",
      source_name        = "CYP2C9 (Missing category in Table 3)"
    )
  )

  covariatesDataExcluded <- list(
    BSA = list(
      description        = "Body surface area (Mosteller formula)",
      units              = "m^2",
      type               = "continuous",
      notes              = "Screened in the Lane 2011 univariate analysis as significant for CL but eliminated in the multivariate / backward-stepwise covariate selection in favour of WGT (Methods 'Covariate selection and models'; Results 'S-Warfarin models')."
    ),
    HEIGHT = list(
      description        = "Body height",
      units              = "m",
      type               = "continuous",
      notes              = "Screened in the Lane 2011 covariate selection but not retained in the final S-warfarin model."
    ),
    CONMED_AMIO = list(
      description        = "Concomitant amiodarone use indicator",
      units              = "(binary)",
      type               = "binary",
      notes              = "Lane 2011 explicitly tested amiodarone (the most commonly used P450 inhibitor in the cohort; 20 of 306 / 6.5% on amiodarone per Table 1) and found no significant impact on the OFV in either univariate or stepwise analysis (Results 'S-Warfarin models'). Notable as a negative finding because amiodarone is widely thought to potentiate warfarin via CYP2C9 inhibition (Discussion paragraph 2)."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 306L,
    n_studies      = 1L,
    age_range      = "19-95 years",
    age_mean       = "66.4 years (Lane 2011 Table 1, total cohort)",
    age_reference  = "69.8 years (model equation reference; not the cohort mean)",
    weight_range   = "36-172 kg",
    weight_mean    = "80.7 kg (Lane 2011 Table 1, total cohort)",
    weight_reference = "70 kg (model equation reference)",
    sex_female_pct = 42,
    race_ethnicity = "UK (Liverpool); race not stratified in the source",
    disease_state  = "Adults on long-term warfarin anticoagulant therapy, indication-unrestricted (Methods 'Patients'): warfarin started for any clinical indication, observational follow-up over 26 weeks",
    dose_range     = "Loading dose (over 3 days) followed by daily maintenance dose individually titrated per UK NHS in-house guidelines; absolute mg-per-day range not reported",
    regions        = "United Kingdom (Royal Liverpool & Broadgreen University Hospital NHS Trust and University Hospital Aintree, Liverpool, 2004-2006)",
    cyp2c9_freq    = "CYP2C9*1/*1 63.7%, *1/*2 19.3%, *1/*3 9.5%, *2/*2 0.3%, *2/*3 2.0%, *3/*3 0.6%, missing 4.6% (Lane 2011 Table 1, S-warfarin cohort n=306)",
    n_pk_records   = "739 S-warfarin plasma concentrations (Lane 2011 Table 1)",
    sampling       = "Sparse sampling at 1, 8, and 26 weeks post-warfarin initiation; samples drawn ~16 h after the patient's previous dose; assay LLOQ 100 ng/mL; assay range 100-5000 ng/mL (Methods 'Determination of plasma warfarin enantiomer concentrations')",
    notes          = "354 patients enrolled (irrespective of indication); 306 included in the S-warfarin PK model and 309 in the R-warfarin model. The two enantiomer models use overlapping but not identical patient subsets and report separate IIV and residual-error magnitudes; modellib('Lane_2011_warfarin_r') is the companion R-warfarin entry."
  )

  ini({
    # ============================================================
    # S-warfarin 1-cmt PK -- Lane 2011 Table 3 (final covariate model)
    # ============================================================
    lcl <- log(0.144) ; label("Apparent oral CL/F (L/h) -- typical, 70 kg woman aged 69.8 y, CYP2C9 *1/*1") # Lane 2011 Table 3 (estimated, SE 0.00647, 95% CI 0.131-0.157)
    lvc <- log(16.6)  ; label("Apparent volume of distribution V/F (L)")                                    # Lane 2011 Table 3 (estimated, SE 1.57,    95% CI 13.5-19.7)
    lka <- fixed(log(1.66)); label("Absorption rate constant Ka (1/h)")                                    # Lane 2011 Table 3 (fixed; sensitivity tested across 1-5 1/h; Methods 'Base models', value from prior literature ref [23])

    # ---- Covariate effects on CL ----
    e_wt_cl  <- 0.321    ; label("Power exponent of body weight on CL (reference 70 kg)")                 # Lane 2011 Table 3 (estimated, SE 0.145, 95% CI 0.037-0.605)
    e_age_cl <- -0.00816 ; label("Fractional linear effect of age on CL (reference 69.8 y, per year)")    # Lane 2011 Table 3 (estimated, SE 0.00206, 95% CI -0.0122 to -0.00412)
    e_male_cl <- 1.12    ; label("CL multiplier for males vs the female reference (unitless)")             # Lane 2011 Table 3 (estimated, SE 0.0705, 95% CI 0.982-1.26)

    # CYP2C9 diplotype CL multipliers (the *1/*1 wild-type reference carries no parameter; multiplier = 1)
    e_cyp2c9_12_cl      <- 0.855 ; label("CL multiplier for CYP2C9 *1/*2 vs *1/*1 reference")            # Lane 2011 Table 3 (estimated, SE 0.0629, 95% CI 0.732-0.978)
    e_cyp2c9_22_cl      <- 0.672 ; label("CL multiplier for CYP2C9 *2/*2 vs *1/*1 reference")            # Lane 2011 Table 3 (estimated, SE 0.0568, 95% CI 0.561-0.783)
    e_cyp2c9_13_cl      <- 0.454 ; label("CL multiplier for CYP2C9 *1/*3 vs *1/*1 reference")            # Lane 2011 Table 3 (estimated, SE 0.0505, 95% CI 0.355-0.553)
    e_cyp2c9_23_cl      <- 0.496 ; label("CL multiplier for CYP2C9 *2/*3 vs *1/*1 reference")            # Lane 2011 Table 3 (estimated, SE 0.121,  95% CI 0.259-0.733)
    e_cyp2c9_33_cl      <- 0.286 ; label("CL multiplier for CYP2C9 *3/*3 vs *1/*1 reference")            # Lane 2011 Table 3 (estimated, SE 0.0324, 95% CI 0.222-0.350)
    e_cyp2c9_missing_cl <- 0.782 ; label("CL multiplier for missing CYP2C9 genotype vs *1/*1 reference")  # Lane 2011 Table 3 (estimated, SE 0.0899, 95% CI 0.606-0.958)

    # ============================================================
    # IIV: block correlation between random effects on log(CL) and log(V)
    #   omega^2 on log scale: omega2 = log(1 + CV^2)
    #     CV(CL) = 41.8% -> omega2_cl = log(1 + 0.418^2) = 0.16113
    #     CV(V)  = 35.8% -> omega2_v  = log(1 + 0.358^2) = 0.12054
    #   Correlation(CL,V) = 0.422 (Lane 2011 Table 3 footnote: "Covariance is expressed as a correlation coefficient")
    #     cov_cl_v = 0.422 * sqrt(0.16113 * 0.12054) = 0.058822
    # ============================================================
    etalcl + etalvc ~ c(0.16113,
                        0.058822, 0.12054)                                                                # Lane 2011 Table 3 (IIV CL 41.8%, IIV V 35.8%, correlation 0.422)

    # ============================================================
    # Residual error
    #   Proportional residual: 31.6% (Lane 2011 Table 3, expressed as approximate CV per footnote)
    #   Additive residual: fixed at 1 mcg/L = 1 ng/mL = 0.001 mg/L (Lane 2011 Table 3 "Additive error 1 Fixed";
    #     concentration units in this model are mg/L per the units list above, so the additive constant is 0.001)
    # ============================================================
    propSd <- 0.316           ; label("Proportional residual SD (fraction)")                              # Lane 2011 Table 3 (estimated, 95% CI 28.5-34.5%)
    addSd  <- fixed(0.001)    ; label("Additive residual SD (mg/L; = 1 ng/mL fixed in the paper)")        # Lane 2011 Table 3 (fixed at low value; the base S-warfarin model estimated 45.8 ng/mL but the final covariate model dropped it to a low fixed value)
  })

  model({
    # ---- Diplotype indicators built from CYP2C9 per-allele counts ----
    not_missing  <- 1 - CYP2C9_MISSING
    is_cyp2c9_11 <- (CYP2C9_S1_COUNT == 2) * not_missing
    is_cyp2c9_12 <- (CYP2C9_S1_COUNT == 1) * (CYP2C9_S2_COUNT == 1) * (CYP2C9_S3_COUNT == 0) * not_missing
    is_cyp2c9_13 <- (CYP2C9_S1_COUNT == 1) * (CYP2C9_S2_COUNT == 0) * (CYP2C9_S3_COUNT == 1) * not_missing
    is_cyp2c9_22 <- (CYP2C9_S2_COUNT == 2) * not_missing
    is_cyp2c9_23 <- (CYP2C9_S1_COUNT == 0) * (CYP2C9_S2_COUNT == 1) * (CYP2C9_S3_COUNT == 1) * not_missing
    is_cyp2c9_33 <- (CYP2C9_S3_COUNT == 2) * not_missing
    is_cyp2c9_missing <- CYP2C9_MISSING

    cyp2c9_mult <- is_cyp2c9_11 +
                   e_cyp2c9_12_cl      * is_cyp2c9_12 +
                   e_cyp2c9_13_cl      * is_cyp2c9_13 +
                   e_cyp2c9_22_cl      * is_cyp2c9_22 +
                   e_cyp2c9_23_cl      * is_cyp2c9_23 +
                   e_cyp2c9_33_cl      * is_cyp2c9_33 +
                   e_cyp2c9_missing_cl * is_cyp2c9_missing

    # ---- Sex multiplier (female reference = 1, male = e_male_cl) ----
    gender_mult <- SEXF + e_male_cl * (1 - SEXF)

    # ---- Individual structural parameters ----
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl *
          (1 + e_age_cl * (AGE - 69.8)) *
          cyp2c9_mult * gender_mult
    vc <- exp(lvc + etalvc)
    ka <- exp(lka)

    # ---- Micro-constants ----
    kel <- cl / vc

    # ---- S-warfarin 1-compartment PK ----
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # ---- Observation and error ----
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
