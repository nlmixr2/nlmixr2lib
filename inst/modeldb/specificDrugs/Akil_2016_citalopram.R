Akil_2016_citalopram <- function() {
  description <- paste(
    "Joint enantiomer-resolved parent-plus-metabolite population PK model for racemic",
    "citalopram in elderly Alzheimer's disease patients treated for agitation in the CitAD",
    "trial (Akil 2016). Four disposition compartments, one per measured analyte -- R-citalopram,",
    "S-citalopram, R-desmethylcitalopram and S-desmethylcitalopram -- fed by a single oral depot",
    "whose first-order absorption rate constant Ka is shared by the two enantiomers and fixed at",
    "1 /h. The racemic capsule splits 50/50 into the two parent compartments. Parent-to-metabolite",
    "conversion is complete, so the apparent metabolic clearance of each parent is also its",
    "formation clearance into the corresponding desmethyl metabolite, and each metabolite is",
    "assumed to share the apparent volume of distribution of its parent enantiomer. Clearance of",
    "the R-enantiomer is slower than that of the S-enantiomer. Covariate effects differ by",
    "enantiomer: R-citalopram apparent metabolic clearance falls with age and is about 30 percent",
    "lower in women; S-citalopram apparent metabolic clearance falls with age, rises with body",
    "weight and is about 36 percent higher in CYP2C19 extensive/rapid metabolizers than in",
    "intermediate/poor metabolizers; both desmethylcitalopram clearances rise with body weight.",
    "All clearances and volumes are apparent (divided by the unknown oral bioavailability F).",
    sep = " "
  )
  reference <- paste(
    "Akil A, Bies RR, Pollock BG, Avramopoulos D, Devanand DP, Mintzer JE, Porsteinsson AP,",
    "Schneider LS, Weintraub D, Yesavage J, Shade DM, Lyketsos CG.",
    "A population pharmacokinetic model for R- and S-citalopram and desmethylcitalopram in",
    "Alzheimer's disease patients with agitation.",
    "J Pharmacokinet Pharmacodyn. 2016 Feb;43(1):99-109.",
    "doi:10.1007/s10928-015-9457-6",
    sep = " "
  )
  vignette <- "Akil_2016_citalopram"

  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ng/mL"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. `depot` holds the racemic capsule contents; the 50/50
  # split into the two parent compartments happens on transfer.
  compartmentData <- list(
    depot = list(
      analyte = "racemic citalopram",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central_r_enant = list(
      analyte = "R-citalopram",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    central_s_enant = list(
      analyte = "S-citalopram",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    central_dcit_r_enant = list(
      analyte = "R-desmethylcitalopram",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    central_dcit_s_enant = list(
      analyte = "S-desmethylcitalopram",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  covariateData <- list(
    AGE = list(
      description = "Subject age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Centred on 60 years, not on the cohort mean: Akil 2016 Methods 'Final model' states",
        "'Age was centered on a value of 60 years', and the Results give the retained",
        "relationships explicitly as CLRp/F = CL0/F x (Age/60)^-0.822 and",
        "CLSp/F = CL0/F x (Age/60)^-1.33. The centring value therefore sits well below the",
        "cohort mean of 77.8 years, so the tabulated CL0 estimates (13 / 9.05 L/h for",
        "R-citalopram, 22.1 / 16.3 / 16.8 L/h for S-citalopram) are extrapolated typical values",
        "at age 60 rather than typical values at the median subject. Evaluating the model at the",
        "cohort mean age reproduces the paper's own post-hoc empirical-Bayes means closely",
        "(10.50 vs 10.59 L/h in men, 7.31 vs 7.25 L/h in women; Akil 2016 Discussion), which is",
        "what confirms the centring. Age was not retained on either desmethylcitalopram",
        "clearance.",
        sep = " "
      ),
      source_name = "Age (Akil 2016 Table 1; 'Age/60' in the Results covariate equations)"
    ),
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Centred on 70 kg (Akil 2016 Methods 'Final model': 'weight was centered on a value of",
        "70 kg'), which is close to the cohort mean of 71.5 kg. Retained on S-citalopram",
        "apparent metabolic clearance as (WT/70)^0.75 and on both desmethylcitalopram apparent",
        "clearances as (WT/70)^0.75 (Akil 2016 Results). Not retained on R-citalopram apparent",
        "metabolic clearance ('No significant effects of body weight, BMI and CYP2C19 genotype",
        "were found on the apparent metabolic clearance of R-citalopram') and not retained on",
        "either apparent volume of distribution. The paper reports the exponent as 0.75 for all",
        "three weight effects and describes the relationship as a 'centered power function' in",
        "the Results text while the Fig. 6c / Fig. 7 captions call the fitted line 'linear'; the",
        "power form in the Results equations is the one implemented (see vignette Errata).",
        sep = " "
      ),
      source_name = "Weight (Akil 2016 Table 1; 'WT/70' in the Results covariate equations)"
    ),
    SEXF = list(
      description = "Female sex indicator; 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "Akil 2016 Methods 'Final model' fits sex as two separate typical values",
        "(theta_1p for females, theta_2p for males) rather than as a reference plus a shift, and",
        "Table 3 reports both: CL Rp /F = 13 L/h for male and 9.05 L/h for female. The model file",
        "encodes male as the reference (SEXF = 0) and carries the female/male ratio 9.05/13 =",
        "0.6962 as e_sexf_cl_r_enant so that both printed estimates remain literally on the",
        "parameter line. Retained on R-citalopram apparent metabolic clearance only: 'No effect",
        "of sex was found on S-citalopram or R,S-desmethylcitalopram apparent clearance'",
        "(Akil 2016 Discussion). The cohort was 40/81 female (49.4 percent).",
        sep = " "
      ),
      source_name = "Sex, Male / Female (Akil 2016 Table 1)"
    ),
    CYP2C19_IM = list(
      description = "CYP2C19 intermediate-metabolizer phenotype indicator; 1 = IM, 0 = otherwise",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (extensive or rapid metabolizer, when CYP2C19_PM and CYP2C19_MISSING are also 0)",
      notes = paste(
        "Akil 2016 Methods 'Final model' regrouped CYP2C19 genotype into three levels",
        "(EM/RM = 1, IM/PM = 2, missing = 3) and estimated one typical clearance per level, so",
        "the intermediate and poor metabolizers share a single pooled coefficient. The general-",
        "scope canonical phenotype columns CYP2C19_IM and CYP2C19_PM are used rather than the",
        "paper-specific composite CYP2C19_NON_EM because (a) they preserve the IM/PM distinction",
        "in the data even though this fit pools it, and (b) CYP2C19_NON_EM's registered reference",
        "category is the homozygous *1/*1 extensive metabolizer, whereas Akil 2016 pools the",
        "rapid metabolizers in with the extensive metabolizers as the reference. Accordingly",
        "e_cyp2c19_im_cl_s_enant and e_cyp2c19_pm_cl_s_enant carry the SAME value -- they are one",
        "estimate, not two. Cohort: 17 IM (21 percent) and 3 PM (3.7 percent) of 81",
        "(Akil 2016 Table 1). Retained on S-citalopram apparent metabolic clearance only; 'No",
        "impact of CYP2C19 genotype was observed in our analysis on the apparent metabolic",
        "clearance of R-citalopram' (Akil 2016 Discussion).",
        sep = " "
      ),
      source_name = "Intermediate metabolizers (Akil 2016 Table 1); 'IM/PM = 2' (Methods 'Final model')"
    ),
    CYP2C19_PM = list(
      description = "CYP2C19 poor-metabolizer phenotype indicator; 1 = PM, 0 = otherwise",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (extensive or rapid metabolizer, when CYP2C19_IM and CYP2C19_MISSING are also 0)",
      notes = paste(
        "Companion to CYP2C19_IM; see that entry for why the two indicators carry a single shared",
        "estimate in this model. Cohort: 3 PM of 81 (3.7 percent; Akil 2016 Table 1).",
        sep = " "
      ),
      source_name = "Poor metabolizers (Akil 2016 Table 1); 'IM/PM = 2' (Methods 'Final model')"
    ),
    CYP2C19_MISSING = list(
      description = "CYP2C19 genotype-missing indicator; 1 = genotype not determined, 0 = genotype known",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (genotype known)",
      notes = paste(
        "Akil 2016 treated an undetermined CYP2C19 genotype as its own covariate level with its",
        "own typical clearance (Methods 'Final model': 'EM/RM = 1, IM/PM = 2 and missing = 3'),",
        "rather than imputing those subjects into the extensive-metabolizer reference. 15 of 81",
        "subjects (18.5 percent) were missing genotype (Akil 2016 Table 1). The fitted",
        "missing-group clearance, 16.8 L/h, lands close to the IM/PM value of 16.3 L/h rather than",
        "to the 22.1 L/h EM/RM reference. This coefficient describes a mixture of unmeasured",
        "phenotypes weighted by their prevalence in the CitAD cohort and does not transfer to",
        "populations with different CYP2C19 allele frequencies. When CYP2C19_MISSING = 1, both",
        "CYP2C19_IM and CYP2C19_PM must be 0.",
        sep = " "
      ),
      source_name = "Missing (Akil 2016 Table 1); 'CL Sp /F for Missing' (Table 3)"
    )
  )

  covariatesDataExcluded <- list(
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Screened as a continuous covariate on every pharmacokinetic parameter using the centred",
        "additive and power models of Akil 2016 Methods 'Final model' ('Both continuous covariates",
        "(age, weight, and BMI) and discrete covariates (CYP2C19 genotype, and sex) were tested'),",
        "centred on 25, and not retained in the final model -- it appears in no row of the",
        "covariate-selection table (Akil 2016 Table 2) and is named explicitly among the",
        "non-significant effects on R-citalopram clearance in the Results. Cohort mean (SD) 26.3",
        "(5.2), range 15.4-41.6. Akil 2016 labels the BMI units 'lbs/in^2' in Table 1 and in",
        "Methods 'Final model'; the reported values are in the conventional kg/m^2 range, so that",
        "label is a typographical error in the source (see vignette Errata).",
        sep = " "
      ),
      source_name = "Body mass index (Akil 2016 Table 1)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 81L,
    n_studies = 1L,
    n_observations = paste(
      "205 R-citalopram, 205 S-citalopram, 179 R-desmethylcitalopram and 109",
      "S-desmethylcitalopram plasma concentrations (Akil 2016 Table 1); 2.5, 2.5, 2.2 and 1.3",
      "observations per subject respectively",
      sep = " "
    ),
    age_range = "mean (SD) 77.8 (8.2) years, range 47-90",
    weight_range = "mean (SD) 71.5 (17.2) kg, range 40-122.3",
    bmi_range = "mean (SD) 26.3 (5.2), range 15.4-41.6",
    sex_female_pct = 49.4,
    race_ethnicity = "Not reported in Akil 2016",
    disease_state = "Alzheimer's disease with clinically significant agitation",
    dose_range = paste(
      "Oral racemic citalopram started at 10 mg once daily and titrated over 2 weeks to a target",
      "of 30 mg once daily, given as three 10 mg capsules in the morning",
      sep = " "
    ),
    genotype = paste(
      "CYP2C19: 43 extensive (53.1 percent), 3 rapid (3.7 percent), 17 intermediate (21 percent),",
      "3 poor (3.7 percent) metabolizers; 15 missing (18.5 percent)",
      sep = " "
    ),
    regions = "United States (CitAD, a multi-centre randomised trial)",
    notes = paste(
      "Sparse sampling: plasma samples were drawn at weeks 3, 6 and 9 of treatment, i.e. at",
      "steady state, with no intensive profile. 94 patients received citalopram and provided",
      "concentration samples; 81 contributed to the population analysis. Concentrations were",
      "measured by chiral HPLC-UV with a limit of quantitation of 5 ng/mL for each analyte except",
      "S-desmethylcitalopram, where it was 10 ng/mL. Because only oral racemate was given, no",
      "bioavailability could be estimated and all clearances and volumes are apparent; Akil 2016",
      "Discussion notes that some of the observed R-versus-S exposure difference may in principle",
      "be a bioavailability rather than a clearance difference. The final model assumes complete",
      "parent-to-metabolite conversion even though 8-10 percent of S-citalopram is excreted",
      "unchanged in urine (Akil 2016 Discussion).",
      sep = " "
    )
  )

  ini({
    # ========================================================================
    # Absorption. Akil 2016 Results: 'The oral absorption rate constant was
    # assumed to be equal for both enantiomers', and Table 3 reports Ka as
    # 1 (Fixed) in both the R- and the S-enantiomer blocks, with bootstrap
    # 'NA'. One shared, fixed Ka therefore serves the whole model.
    # ========================================================================
    lka <- fixed(log(1)); label("Absorption rate constant Ka, shared by both enantiomers (1/h)")

    # ========================================================================
    # R-enantiomer disposition. Akil 2016 Table 3, 'R-enantiomer' block,
    # 'Final model estimate' column. CL Rp /F is the apparent METABOLIC
    # clearance of R-citalopram: under the final model's complete
    # parent-to-metabolite conversion it is simultaneously the elimination
    # clearance of R-citalopram and the formation clearance of
    # R-desmethylcitalopram (Akil 2016 equations (1) and (3)).
    # The two sex-specific estimates are encoded as a male reference plus the
    # printed female/male ratio; see e_sexf_cl_r_enant below.
    # ========================================================================
    lcl_r_enant      <- log(13);   label("Apparent metabolic clearance of R-citalopram CLRp/F in men at age 60 (L/h)")   # Akil 2016 Table 3: 'CL Rp /F for male, L/h' = 13 (bootstrap median 13.8, 95% CI 13.5-14.1)
    lvc_r_enant      <- log(1830); label("Apparent volume of distribution of the R-enantiomer V/F (L)")                  # Akil 2016 Table 3: 'V/F, L' = 1830 (bootstrap median 1605, 95% CI 1440-2090); shared by R-citalopram and R-desmethylcitalopram
    lcl_dcit_r_enant <- log(24.4); label("Apparent clearance of R-desmethylcitalopram CLRm/F at 70 kg (L/h)")            # Akil 2016 Table 3: 'CL Rm /F, L/h' = 24.4 (bootstrap median 23.5, 95% CI 23.1-23.7)

    # ========================================================================
    # S-enantiomer disposition. Akil 2016 Table 3, 'S-enantiomer' block.
    # The three CYP2C19 strata are encoded as an EM/RM reference plus the
    # printed stratum ratios; see the CYP2C19 effects below.
    # ========================================================================
    lcl_s_enant      <- log(22.1); label("Apparent metabolic clearance of S-citalopram CLSp/F in CYP2C19 EM/RM at age 60 and 70 kg (L/h)")  # Akil 2016 Table 3: 'CL Sp /F for EM/RM, L/h' = 22.1 (bootstrap median 21.9, 95% CI 21.2-22.8)
    lvc_s_enant      <- log(1390); label("Apparent volume of distribution of the S-enantiomer V/F (L)")                                     # Akil 2016 Table 3: 'V/F, L' = 1390 (bootstrap median 1310, 95% CI 1130-1420); shared by S-citalopram and S-desmethylcitalopram
    lcl_dcit_s_enant <- log(38.8); label("Apparent clearance of S-desmethylcitalopram CLSm/F at 70 kg (L/h)")                               # Akil 2016 Table 3: 'CL Sm /F, L/h' = 38.8 (bootstrap median 38.9, 95% CI 38.4-39.2)

    # ========================================================================
    # Covariate effects.
    #
    # Continuous covariates use the centred power form Akil 2016 selected in
    # the Results ('A centered power model was chosen to model the effects of
    # continuous covariates (age and weight) on pharmacokinetic parameter
    # estimates'), with the centring constants given in Methods 'Final model'
    # (age 60 years, weight 70 kg).
    #
    # Categorical covariates were fitted by Akil 2016 as one typical value per
    # level (Methods 'Final model': theta_1p for females and theta_2p for
    # males; theta_1p for EM/RM, theta_2p for IM/PM and theta_3p for missing).
    # They are re-expressed here as a reference level plus the ratio of the two
    # printed estimates, so that every number Akil 2016 published stays
    # literally on the parameter line and the reference category is explicit.
    # ========================================================================
    e_age_cl_r_enant <- -0.822; label("Power exponent on (AGE/60) for R-citalopram CLRp/F (unitless)")   # Akil 2016 Results: 'CLRp/F = CL0/F x (Age/60)^-0.822'
    e_age_cl_s_enant <- -1.33;  label("Power exponent on (AGE/60) for S-citalopram CLSp/F (unitless)")   # Akil 2016 Results: 'CLSp/F = CL0/F x (Age/60)^-1.33'
    e_wt_cl_s_enant  <- 0.75;   label("Power exponent on (WT/70) for S-citalopram CLSp/F (unitless)")    # Akil 2016 Results: 'CLSp/F = CL0/F x (WT/70)^0.75'

    e_wt_cl_dcit_r_enant <- 0.75; label("Power exponent on (WT/70) for R-desmethylcitalopram CLRm/F (unitless)")  # Akil 2016 Results: 'CLm/F = CL0/F x (WT/70)^0.75', stated once for both desmethylcitalopram enantiomers
    e_wt_cl_dcit_s_enant <- 0.75; label("Power exponent on (WT/70) for S-desmethylcitalopram CLSm/F (unitless)")  # Akil 2016 Results: 'CLm/F = CL0/F x (WT/70)^0.75', stated once for both desmethylcitalopram enantiomers

    e_sexf_cl_r_enant <- 9.05 / 13; label("Female/male ratio of R-citalopram CLRp/F (unitless)")  # Akil 2016 Table 3: 9.05 L/h female / 13 L/h male = 0.6962; Results quote the same contrast as 'approximately 30 % higher in males'

    # Akil 2016 pooled the intermediate and poor metabolizers into one stratum
    # (Methods 'Final model': 'IM/PM = 2'), so the next two lines carry ONE
    # estimate applied through two canonical phenotype indicators -- they are
    # deliberately equal, not two independent numbers.
    e_cyp2c19_im_cl_s_enant      <- 16.3 / 22.1; label("IM/EM-RM ratio of S-citalopram CLSp/F (unitless)")       # Akil 2016 Table 3: 16.3 L/h IM/PM / 22.1 L/h EM/RM = 0.7376; Results quote the reciprocal contrast as 'about 36 % higher'
    e_cyp2c19_pm_cl_s_enant      <- 16.3 / 22.1; label("PM/EM-RM ratio of S-citalopram CLSp/F (unitless)")       # Akil 2016 Table 3: same pooled IM/PM estimate as e_cyp2c19_im_cl_s_enant
    e_cyp2c19_missing_cl_s_enant <- 16.8 / 22.1; label("Missing-genotype/EM-RM ratio of S-citalopram CLSp/F (unitless)")  # Akil 2016 Table 3: 'CL Sp /F for Missing, L/h' = 16.8 / 22.1 = 0.7602 (the Results text prints 16.6; see vignette Errata)

    # ========================================================================
    # Between-subject variability. Akil 2016 Methods 'Base model' specifies a
    # log-normal BSV model, Pj = PTV x exp(eta_p). Table 3 reports each omega
    # in PERCENT, so the internal variance is omega^2 = log(1 + CV^2). The
    # column's 'Variance of the BSV of ...' description text is a mislabel: a
    # variance cannot carry a percent unit, and the sibling residual-error row
    # in the same column is labelled 'ng/ml' -- an SD unit, not the (ng/ml)^2
    # a variance would need. Reading the column as CV% also makes the
    # covariance row admissible (see below) and reproduces the paper's own
    # empirical-Bayes spreads; reading it as a variance does neither.
    # See the vignette Errata for the full argument.
    # ========================================================================
    etalcl_r_enant      ~ 0.0672758  # Akil 2016 Table 3, R block, 'x CLp , %' = 26.38 (bootstrap median 28.7, 95% CI 27.8-30); log(1 + 0.2638^2) = 0.0672758
    etalvc_r_enant      ~ 1.3296947  # Akil 2016 Table 3, R block, 'x V , %' = 166.73 (bootstrap median 107.2, 95% CI 92.6-121.2); log(1 + 1.6673^2) = 1.3296947
    etalcl_dcit_r_enant ~ 0.0895639  # Akil 2016 Table 3, R block, 'x CLm , %' = 30.61 (bootstrap median 34.9, 95% CI 34.5-35.9); log(1 + 0.3061^2) = 0.0895639

    # Akil 2016 Table 3 reports a covariance between the S-enantiomer
    # clearance and volume etas, 'x CLp ; V , %' = 47.1. On the CV% reading of
    # this column the entry can only be a CORRELATION: taken as a raw
    # covariance it implies corr = 0.471 / (0.3703 x 0.6707) = 1.90, which is
    # outside [-1, 1] and would make the omega block non-positive-definite.
    # cov = 0.471 x sqrt(0.1371460) x sqrt(0.4498415) = 0.1169882.
    etalcl_s_enant + etalvc_s_enant ~ c(0.1371460,
                                        0.1169882, 0.4498415)  # Akil 2016 Table 3, S block: 'x CLp , %' = 38.34, 'x CLp ; V , %' = 47.1, 'x V , %' = 75.37
    etalcl_dcit_s_enant ~ 0.0411266  # Akil 2016 Table 3, S block, 'x CLm , %' = 20.49 (bootstrap median 20.1, 95% CI 19.4-20.6); log(1 + 0.2049^2) = 0.0411266

    # ========================================================================
    # Residual error. Akil 2016 Results: 'The residual error was separated with
    # an additive structure for R-citalopram, a proportional structure for
    # R-desmethylcitalopram and both S-citalopram and desmethylcitalopram.'
    # That sentence accounts for exactly the three sigma rows in Table 3: the
    # R block's additive (ng/ml) row belongs to R-citalopram and its
    # proportional (%) row to R-desmethylcitalopram, while the S block's single
    # proportional row is shared by S-citalopram and S-desmethylcitalopram.
    # ========================================================================
    addSd_r_enant      <- 13.42;  label("Additive residual SD on R-citalopram (ng/mL)")                          # Akil 2016 Table 3, R block: 'r , ng/ml (additive)' = 13.42 (bootstrap median 13.6, 95% CI 13.3-13.9)
    propSd_dcit_r_enant <- 0.2154; label("Proportional residual SD on R-desmethylcitalopram (fraction)")         # Akil 2016 Table 3, R block: 'r , % (proportional)' = 21.54 (bootstrap median 20.7, 95% CI 20.4-21.3)
    propSd_s_enant     <- 0.2161; label("Proportional residual SD on S-citalopram (fraction)")                   # Akil 2016 Table 3, S block: 'r , % (proportional)' = 21.61 (bootstrap median 21.6, 95% CI 21-22.1)
    propSd_dcit_s_enant <- 0.2161; label("Proportional residual SD on S-desmethylcitalopram (fraction)")         # Akil 2016 Table 3, S block: the single S proportional sigma, 21.61, covers 'both S-citalopram and desmethylcitalopram' (Results); repeated here because rxode2 needs one named SD per endpoint
  })

  model({
    # ---- Absorption --------------------------------------------------------
    ka <- exp(lka)

    # ---- R-enantiomer individual parameters --------------------------------
    # Akil 2016 Results: CLRp/F = CL0/F x (Age/60)^-0.822, with CL0/F differing
    # by sex (13 L/h male, 9.05 L/h female).
    cl_r_enant <- exp(lcl_r_enant + etalcl_r_enant) *
      e_sexf_cl_r_enant^SEXF *
      (AGE / 60)^e_age_cl_r_enant
    vc_r_enant <- exp(lvc_r_enant + etalvc_r_enant)
    cl_dcit_r_enant <- exp(lcl_dcit_r_enant + etalcl_dcit_r_enant) *
      (WT / 70)^e_wt_cl_dcit_r_enant

    # ---- S-enantiomer individual parameters --------------------------------
    # The three CYP2C19 levels are mutually exclusive, so this evaluates to 1
    # in the EM/RM reference group and to the corresponding stratum ratio
    # otherwise. Akil 2016 Results:
    #   CLSp/F = CL0/F x (Age/60)^-1.33 x (WT/70)^0.75,
    # with CL0/F = 22.1 (EM/RM), 16.3 (IM/PM) or 16.8 (missing) L/h.
    cyp2c19_cl_s_enant <- 1 +
      (e_cyp2c19_im_cl_s_enant - 1) * CYP2C19_IM +
      (e_cyp2c19_pm_cl_s_enant - 1) * CYP2C19_PM +
      (e_cyp2c19_missing_cl_s_enant - 1) * CYP2C19_MISSING

    cl_s_enant <- exp(lcl_s_enant + etalcl_s_enant) *
      cyp2c19_cl_s_enant *
      (AGE / 60)^e_age_cl_s_enant *
      (WT / 70)^e_wt_cl_s_enant
    vc_s_enant <- exp(lvc_s_enant + etalvc_s_enant)
    cl_dcit_s_enant <- exp(lcl_dcit_s_enant + etalcl_dcit_s_enant) *
      (WT / 70)^e_wt_cl_dcit_s_enant

    # ---- Micro-constants ---------------------------------------------------
    # Akil 2016 assumed 'Both the parent and metabolite for the two enantiomers
    # ... have the same volume of distribution', which is why the metabolite
    # elimination rate constants below divide by the PARENT enantiomer's
    # volume: equations (3) and (4) write the metabolite loss term as
    # (CLRm/VR) x C(3) and (CLSm/VS) x C(4), not CLRm/VRm.
    kel_r_enant      <- cl_r_enant      / vc_r_enant
    kel_s_enant      <- cl_s_enant      / vc_s_enant
    kel_dcit_r_enant <- cl_dcit_r_enant / vc_r_enant
    kel_dcit_s_enant <- cl_dcit_s_enant / vc_s_enant

    # ---- Disposition -------------------------------------------------------
    # Akil 2016 equations (1)-(4), written by the authors in concentration
    # units with the absorption input abbreviated 'Ka x Dose / V'; recast here
    # in amounts (mg), which is the equivalent form rxode2 solves and which
    # makes the complete parent-to-metabolite conversion mass-conserving. The
    # authors' C(1)..C(4) map to Cc_r_enant, Cc_s_enant, Cc_dcit_r_enant and
    # Cc_dcit_s_enant respectively.
    #
    # Citalopram is given as a single capsule of the 50/50 racemate ('both
    # enantiomers are administered simultaneously in one pill containing the
    # racemic citalopram' ... 'citalopram is administered as 50/50 racemic
    # mixture', Akil 2016 Discussion), so the dose record carries the RACEMIC
    # amount -- 30 mg for the CitAD target dose -- and half of it reaches each
    # parent compartment. The two enantiomers share one depot because Ka is
    # common to both; a pair of half-dosed depots would be identical.
    #
    # That the tabulated apparent clearances are scaled to the per-enantiomer
    # (half) dose rather than the racemic dose is confirmed non-circularly by
    # the paper's own exposure summary: Dose/CL using the empirical-Bayes mean
    # clearances and a per-enantiomer dose gives 1.5-1.7 mg*h/L for
    # R-citalopram and 0.9-1.1 for S-citalopram against the reported
    # 1.46 +/- 0.58 and 0.97 +/- 0.45 mg*h/L (Akil 2016 Discussion), whereas
    # the full racemic dose would give roughly twice those values.
    d/dt(depot)                <- -ka * depot
    d/dt(central_r_enant)      <-  0.5 * ka * depot - kel_r_enant * central_r_enant
    d/dt(central_s_enant)      <-  0.5 * ka * depot - kel_s_enant * central_s_enant
    d/dt(central_dcit_r_enant) <-  kel_r_enant * central_r_enant -
      kel_dcit_r_enant * central_dcit_r_enant
    d/dt(central_dcit_s_enant) <-  kel_s_enant * central_s_enant -
      kel_dcit_s_enant * central_dcit_s_enant

    # ---- Observations ------------------------------------------------------
    # Amounts are in mg and volumes in L, so amount/volume is mg/L = ug/mL;
    # multiply by 1000 to report ng/mL, the unit in which the limits of
    # quantitation (5 and 10 ng/mL) and the additive residual SD (13.42 ng/mL)
    # are expressed.
    Cc_r_enant      <- 1000 * central_r_enant      / vc_r_enant
    Cc_s_enant      <- 1000 * central_s_enant      / vc_s_enant
    Cc_dcit_r_enant <- 1000 * central_dcit_r_enant / vc_r_enant
    Cc_dcit_s_enant <- 1000 * central_dcit_s_enant / vc_s_enant

    Cc_r_enant      ~ add(addSd_r_enant)
    Cc_s_enant      ~ prop(propSd_s_enant)
    Cc_dcit_r_enant ~ prop(propSd_dcit_r_enant)
    Cc_dcit_s_enant ~ prop(propSd_dcit_s_enant)
  })
}
