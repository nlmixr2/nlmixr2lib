Lane_2011_warfarin_r <- function() {
  description <- "R-warfarin population PK (1-compartment, first-order absorption) in adults on long-term warfarin therapy (Lane 2011). Bodyweight, age, CYP2C19 rs3814637 genotype, and CYP3A4 rs2242480 (CYP3A4*1G) genotype influence apparent clearance; volume of distribution carries no covariates. Block correlation between random effects on CL and V. S-warfarin is reported separately in the same paper (modellib('Lane_2011_warfarin_s'))."
  reference <- paste(
    "Lane S, Al-Zubiedi S, Hatch E, Matthews I, Jorgensen AL, Deloukas P, Daly AK,",
    "Park BK, Aarons L, Ogungbenro K, Kamali F, Hughes D, Pirmohamed M.",
    "The population pharmacokinetics of R- and S-warfarin: effect of genetic and",
    "clinical factors.",
    "Br J Clin Pharmacol. 2012;73(1):66-76.",
    "doi:10.1111/j.1365-2125.2011.04051.x.",
    "PMID: 21692829.",
    "PK parameters and weight / age / CYP2C19 (rs3814637) / CYP3A4 (rs2242480)",
    "effects from Table 3 (final covariate model); structural-equation form",
    "(CL_i = theta_CL * (WT/70)^theta_wgt * (1 + theta_age*(AGE-69.8)) * theta_CYP2C19",
    "* theta_CYP3A4 * exp(eta_CL)) from the Results paragraph following Tables 2 and 3."
  )
  vignette <- "Lane_2011_warfarin"

  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at baseline. Power scaling on CL with reference weight 70 kg (Lane 2011 final R-warfarin equation; the typical CL of 0.125 L/h is reported for a 70-kg, 69.8-year-old subject with CYP2C19 wild-homozygote and CYP3A4 wild-homozygote). Cohort mean (range) body weight 80.7 (36-172) kg (Lane 2011 Table 1).",
      source_name        = "WGT"
    ),
    AGE = list(
      description        = "Age at baseline",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at baseline. Linear fractional effect on CL with reference age 69.8 years (Lane 2011 final R-warfarin equation). Cohort mean (range) age 66.4 (19-95) years (Lane 2011 Table 1); 69.8 y is the equation reference rather than the cohort mean.",
      source_name        = "AGE"
    ),
    SNP_CYP2C19_RS3814637_VAR_COUNT = list(
      description        = "CYP2C19 rs3814637 variant-allele count (0 = wild-homozygote, 1 = heterozygote, 2 = variant homozygote)",
      units              = "(count, 0/1/2)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-invariant germline genotype. Paired with SNP_CYP2C19_RS3814637_MISSING. Cohort distribution (Lane 2011 Table 1, R-warfarin subset n=309): wild-homozygote 70.9%, heterozygote 8.1%, mutant-homozygote 1.3%, missing 19.7%. The rs3814637 SNP is NOT in linkage disequilibrium with the canonical CYP2C19*2 loss-of-function allele (rs4244285) in this patient population per Lane 2011 Discussion paragraph 'R-Warfarin models'.",
      source_name        = "CYP2C19 (rs3814637) -- coded as wild-type / heterozygote / mutant-homozygote in Lane 2011 Methods"
    ),
    SNP_CYP2C19_RS3814637_MISSING = list(
      description        = "Binary indicator: CYP2C19 rs3814637 genotype not measured / not available",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CYP2C19 rs3814637 genotype known)",
      notes              = "Lane 2011 fits a separate CL multiplier (0.804) for the missing-genotype subgroup (61 of 309 R-warfarin subjects). When SNP_CYP2C19_RS3814637_MISSING == 1 the VAR_COUNT column should be 0 so the genotype-indicator products evaluate to zero.",
      source_name        = "CYP2C19 (Missing category in Table 3)"
    ),
    SNP_CYP3A4_RS2242480_VAR_COUNT = list(
      description        = "CYP3A4 rs2242480 variant-allele count (0 = wild-homozygote, 1 = heterozygote, 2 = variant homozygote)",
      units              = "(count, 0/1/2)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-invariant germline genotype. Paired with SNP_CYP3A4_RS2242480_MISSING. Cohort distribution (Lane 2011 Table 1, R-warfarin subset n=309): wild-homozygote 73.1%, heterozygote 13.6%, mutant-homozygote 1.0%, missing 12.3%. rs2242480 contributes to the CYP3A4*1G haplotype; its functional effect remains controversial per Lane 2011 Discussion paragraph 'R-Warfarin models'.",
      source_name        = "CYP3A4 (rs2242480) -- coded as wild-type / heterozygote / mutant-homozygote in Lane 2011 Methods"
    ),
    SNP_CYP3A4_RS2242480_MISSING = list(
      description        = "Binary indicator: CYP3A4 rs2242480 genotype not measured / not available",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CYP3A4 rs2242480 genotype known)",
      notes              = "Lane 2011 fits a separate CL multiplier (0.937) for the missing-genotype subgroup (38 of 309 R-warfarin subjects). When SNP_CYP3A4_RS2242480_MISSING == 1 the VAR_COUNT column should be 0 so the genotype-indicator products evaluate to zero.",
      source_name        = "CYP3A4 (Missing category in Table 3)"
    )
  )

  covariatesDataExcluded <- list(
    BSA = list(
      description        = "Body surface area (Mosteller formula)",
      units              = "m^2",
      type               = "continuous",
      notes              = "Screened in the Lane 2011 R-warfarin univariate analysis but not retained in the final model (Methods 'Covariate selection and models')."
    ),
    HEIGHT = list(
      description        = "Body height",
      units              = "m",
      type               = "continuous",
      notes              = "Screened but not retained in the final R-warfarin model."
    ),
    SEXF = list(
      description        = "Sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      notes              = "Screened but not retained as significant in the R-warfarin model (unlike S-warfarin where sex IS retained). See modellib('Lane_2011_warfarin_s') for the S-warfarin extraction that does include SEXF."
    ),
    CONMED_AMIO = list(
      description        = "Concomitant amiodarone use indicator",
      units              = "(binary)",
      type               = "binary",
      notes              = "Screened but not retained in the R-warfarin model."
    ),
    CYP2C9_S1_COUNT = list(
      description        = "Count of CYP2C9*1 alleles (paired with S2 and S3 counts)",
      units              = "(count, 0/1/2)",
      type               = "continuous",
      notes              = "CYP2C9 genotype was tested in the R-warfarin model but not retained. CYP2C9 IS retained in the S-warfarin model (modellib('Lane_2011_warfarin_s')); only the S-enantiomer is metabolised primarily by CYP2C9 (Lane 2011 Discussion paragraph 'R-Warfarin models'). CYP1A2 SNPs (three rsids; specific identifiers not enumerated in the Lane 2011 paper) were also screened for R-warfarin and found nonsignificant; they are not registered here because the source does not name them."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 309L,
    n_studies      = 1L,
    age_range      = "19-95 years",
    age_mean       = "66.4 years (Lane 2011 Table 1, total cohort)",
    age_reference  = "69.8 years (model equation reference; not the cohort mean)",
    weight_range   = "36-172 kg",
    weight_mean    = "80.7 kg (Lane 2011 Table 1, total cohort)",
    weight_reference = "70 kg (model equation reference)",
    sex_female_pct = 41,
    race_ethnicity = "UK (Liverpool); race not stratified in the source",
    disease_state  = "Adults on long-term warfarin anticoagulant therapy, indication-unrestricted (Methods 'Patients'): warfarin started for any clinical indication, observational follow-up over 26 weeks",
    dose_range     = "Loading dose (over 3 days) followed by daily maintenance dose individually titrated per UK NHS in-house guidelines; absolute mg-per-day range not reported",
    regions        = "United Kingdom (Royal Liverpool & Broadgreen University Hospital NHS Trust and University Hospital Aintree, Liverpool, 2004-2006)",
    cyp2c19_freq   = "rs3814637 wild-homozygote 70.9%, heterozygote 8.1%, mutant-homozygote 1.3%, missing 19.7% (Lane 2011 Table 1, R-warfarin cohort n=309)",
    cyp3a4_freq    = "rs2242480 wild-homozygote 73.1%, heterozygote 13.6%, mutant-homozygote 1.0%, missing 12.3% (Lane 2011 Table 1, R-warfarin cohort n=309)",
    n_pk_records   = "759 R-warfarin plasma concentrations (Lane 2011 Table 1)",
    sampling       = "Sparse sampling at 1, 8, and 26 weeks post-warfarin initiation; samples drawn ~16 h after the patient's previous dose; assay LLOQ 100 ng/mL; assay range 100-5000 ng/mL (Methods 'Determination of plasma warfarin enantiomer concentrations')",
    notes          = "354 patients enrolled (irrespective of indication); 309 included in the R-warfarin PK model and 306 in the S-warfarin model. The two enantiomer models use overlapping but not identical patient subsets and report separate IIV and residual-error magnitudes; modellib('Lane_2011_warfarin_s') is the companion S-warfarin entry."
  )

  ini({
    # ============================================================
    # R-warfarin 1-cmt PK -- Lane 2011 Table 3 (final covariate model)
    # ============================================================
    lcl <- log(0.125) ; label("Apparent oral CL/F (L/h) -- typical, 70 kg subject aged 69.8 y, CYP2C19 wild-hom + CYP3A4 wild-hom") # Lane 2011 Table 3 (estimated, SE 0.00528, 95% CI 0.115-0.135)
    lvc <- log(10.9)  ; label("Apparent volume of distribution V/F (L)")                                    # Lane 2011 Table 3 (estimated, SE 1.16, 95% CI 8.63-13.2)
    lka <- fixed(log(1.66)); label("Absorption rate constant Ka (1/h)")                                    # Lane 2011 Table 3 (fixed; sensitivity tested across 1-5 1/h; Methods 'Base models', value from prior literature ref [23])

    # ---- Covariate effects on CL ----
    e_wt_cl  <- 0.650    ; label("Power exponent of body weight on CL (reference 70 kg)")                 # Lane 2011 Table 3 (estimated, SE 0.132, 95% CI 0.391-0.909)
    e_age_cl <- -0.00657 ; label("Fractional linear effect of age on CL (reference 69.8 y, per year)")    # Lane 2011 Table 3 (estimated, SE 0.00223, 95% CI -0.0109 to -0.00220)

    # CYP2C19 rs3814637 categorical CL multipliers (wild-homozygote reference carries no parameter; multiplier = 1)
    e_cyp2c19_het_cl     <- 0.761 ; label("CL multiplier for CYP2C19 rs3814637 heterozygote vs wild-hom reference")        # Lane 2011 Table 3 (estimated, SE 0.0576, 95% CI 0.648-0.874)
    e_cyp2c19_varhom_cl  <- 0.494 ; label("CL multiplier for CYP2C19 rs3814637 variant homozygote vs wild-hom reference")  # Lane 2011 Table 3 (estimated, SE 0.191, 95% CI 0.120-0.868)
    e_cyp2c19_missing_cl <- 0.804 ; label("CL multiplier for missing CYP2C19 rs3814637 genotype vs wild-hom reference")    # Lane 2011 Table 3 (estimated, SE 0.0523, 95% CI 0.701-0.907)

    # CYP3A4 rs2242480 categorical CL multipliers (wild-homozygote reference carries no parameter; multiplier = 1)
    e_cyp3a4_het_cl      <- 1.32  ; label("CL multiplier for CYP3A4 rs2242480 heterozygote vs wild-hom reference")          # Lane 2011 Table 3 (estimated, SE 0.0978, 95% CI 1.13-1.51)
    e_cyp3a4_varhom_cl   <- 1.06  ; label("CL multiplier for CYP3A4 rs2242480 variant homozygote vs wild-hom reference")    # Lane 2011 Table 3 (estimated, SE 0.172, 95% CI 0.723-1.40)
    e_cyp3a4_missing_cl  <- 0.937 ; label("CL multiplier for missing CYP3A4 rs2242480 genotype vs wild-hom reference")      # Lane 2011 Table 3 (estimated, SE 0.0845, 95% CI 0.771-1.10)

    # ============================================================
    # IIV: block correlation between random effects on log(CL) and log(V)
    #   omega^2 on log scale: omega2 = log(1 + CV^2)
    #     CV(CL) = 43.0% -> omega2_cl = log(1 + 0.430^2) = 0.16975
    #     CV(V)  = 38.3% -> omega2_v  = log(1 + 0.383^2) = 0.13692
    #   Correlation(CL,V) = 0.352 (Lane 2011 Table 3 footnote: "Covariance is expressed as a correlation coefficient")
    #     cov_cl_v = 0.352 * sqrt(0.16975 * 0.13692) = 0.053666
    # ============================================================
    etalcl + etalvc ~ c(0.16975,
                        0.053666, 0.13692)                                                                # Lane 2011 Table 3 (IIV CL 43.0%, IIV V 38.3%, correlation 0.352)

    # ============================================================
    # Residual error
    #   Proportional residual: 31.9% (Lane 2011 Table 3, expressed as approximate CV per footnote)
    #   Additive residual: fixed at 1 ng/mL = 0.001 mg/L (Lane 2011 Table 3 "Additive error 1 Fixed")
    # ============================================================
    propSd <- 0.319           ; label("Proportional residual SD (fraction)")                              # Lane 2011 Table 3 (estimated, 95% CI 29.2-34.5%)
    addSd  <- fixed(0.001)    ; label("Additive residual SD (mg/L; = 1 ng/mL fixed in the paper)")        # Lane 2011 Table 3 (fixed at low value; the additive component dropped out of the final R-warfarin model)
  })

  model({
    # ---- CYP2C19 rs3814637 indicators ----
    not_missing_c19 <- 1 - SNP_CYP2C19_RS3814637_MISSING
    is_cyp2c19_wild    <- (SNP_CYP2C19_RS3814637_VAR_COUNT == 0) * not_missing_c19
    is_cyp2c19_het     <- (SNP_CYP2C19_RS3814637_VAR_COUNT == 1) * not_missing_c19
    is_cyp2c19_varhom  <- (SNP_CYP2C19_RS3814637_VAR_COUNT == 2) * not_missing_c19
    is_cyp2c19_missing <- SNP_CYP2C19_RS3814637_MISSING

    cyp2c19_mult <- is_cyp2c19_wild +
                    e_cyp2c19_het_cl     * is_cyp2c19_het +
                    e_cyp2c19_varhom_cl  * is_cyp2c19_varhom +
                    e_cyp2c19_missing_cl * is_cyp2c19_missing

    # ---- CYP3A4 rs2242480 indicators ----
    not_missing_3a4 <- 1 - SNP_CYP3A4_RS2242480_MISSING
    is_cyp3a4_wild    <- (SNP_CYP3A4_RS2242480_VAR_COUNT == 0) * not_missing_3a4
    is_cyp3a4_het     <- (SNP_CYP3A4_RS2242480_VAR_COUNT == 1) * not_missing_3a4
    is_cyp3a4_varhom  <- (SNP_CYP3A4_RS2242480_VAR_COUNT == 2) * not_missing_3a4
    is_cyp3a4_missing <- SNP_CYP3A4_RS2242480_MISSING

    cyp3a4_mult <- is_cyp3a4_wild +
                   e_cyp3a4_het_cl     * is_cyp3a4_het +
                   e_cyp3a4_varhom_cl  * is_cyp3a4_varhom +
                   e_cyp3a4_missing_cl * is_cyp3a4_missing

    # ---- Individual structural parameters ----
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl *
          (1 + e_age_cl * (AGE - 69.8)) *
          cyp2c19_mult * cyp3a4_mult
    vc <- exp(lvc + etalvc)
    ka <- exp(lka)

    # ---- Micro-constants ----
    kel <- cl / vc

    # ---- R-warfarin 1-compartment PK ----
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # ---- Observation and error ----
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
