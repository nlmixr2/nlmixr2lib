Huang_2022_tacrolimus <- function() {
  description <- paste0(
    "One-compartment population PK model with first-order absorption and ",
    "first-order elimination for oral tacrolimus whole-blood concentrations ",
    "in Chinese children with steroid-resistant or steroid-dependent ",
    "(refractory) nephrotic syndrome (Huang 2022). The absorption rate ",
    "constant ka is fixed at 4.48 1/h from earlier paediatric tacrolimus ",
    "models. Apparent oral clearance CL/F is a power function of age ",
    "(reference 5.3 years) with exponential effects of the CTLA4 rs4553808 ",
    "GA and AA genotypes (reference GG), CYP3A5*3/*3 non-expresser status ",
    "and concomitant Wuzhi capsule. Apparent volume V/F has no covariates. ",
    "Exponential IIV on CL/F and V/F; proportional residual error."
  )
  reference <- paste0(
    "Huang Q, Lin X, Wang Y, Chen X, Zheng W, Zhong X, Shang D, Huang M, ",
    "Gao X, Deng H, Li J, Zeng F, Mo X. Tacrolimus pharmacokinetics in ",
    "pediatric nephrotic syndrome: A combination of population ",
    "pharmacokinetic modelling and machine learning approaches to improve ",
    "individual prediction. Front Pharmacol. 2022;13:942129. ",
    "doi:10.3389/fphar.2022.942129"
  )
  vignette <- "Huang_2022_tacrolimus"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "tacrolimus", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "tacrolimus", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Power effect on CL/F, (AGE/5.3)^0.31 (Eq. 2); 5.3 years is the ",
        "cohort median (Table 1, range 1.1-15.6 years), consistent with the ",
        "Methods covariate form (COV/COVmedian)^theta."
      ),
      source_name = "AGE"
    ),
    SNP_CTLA4_RS4553808_GA = list(
      description = "CTLA4 rs4553808 heterozygous GA genotype indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 together with SNP_CTLA4_RS4553808_AA = 0, i.e. the GG genotype",
      notes = paste0(
        "Text below Eq. 2: 'CTLA4 GA = 1 ... when genotype is CTLA4 GA, 0 ",
        "otherwise'. Enters CL/F as exp(-0.34 x CTLA4 GA). Table 2: GA 29 of ",
        "139 (20.9%); the GG reference group holds only 2 patients (1.4%)."
      ),
      source_name = "CTLA4 GA"
    ),
    SNP_CTLA4_RS4553808_AA = list(
      description = "CTLA4 rs4553808 homozygous AA genotype indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 together with SNP_CTLA4_RS4553808_GA = 0, i.e. the GG genotype",
      notes = paste0(
        "Text below Eq. 2: 'CTLA4 AA = 1 when genotype is ... CTLA4 AA, 0 ",
        "otherwise'. Enters CL/F as exp(-0.15 x CTLA4 AA). Table 2: AA 108 of ",
        "139 (77.7%)."
      ),
      source_name = "CTLA4 AA"
    ),
    CYP3A5_EXPR = list(
      description = "CYP3A5 expresser status (1 = at least one CYP3A5*1 allele)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (CYP3A5*3/*3 nonexpresser)",
      notes = paste0(
        "VALUE INVERSION relative to the source: Huang 2022 codes ",
        "'CYP3A5 *3/*3 = 1 when genotype is CYP3A5 *3/*3, 0 otherwise' ",
        "(text below Eq. 2), so CYP3A5_EXPR = 1 - source indicator and the ",
        "model applies the published coefficient -0.25 to (1 - CYP3A5_EXPR). ",
        "rs776746 by PCR-RFLP. Table 2: *1/*1 18, *1/*3 58, *3/*3 63 of 139."
      ),
      source_name = "CYP3A5 *3/*3"
    ),
    CONMED_WUZHI = list(
      description = "Concomitant Wuzhi capsule indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no Wuzhi capsule)",
      notes = paste0(
        "Text below Eq. 2: 'combWZ = 1 when combined with Wuzhi capsules, 0 ",
        "otherwise'. Enters CL/F as exp(-0.34 x combWZ). Table 1: 6 of 139 ",
        "patients. The paper does not state whether the indicator varied ",
        "within a patient."
      ),
      source_name = "combWZ"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Sex (1 = female)",
      units = "(binary)",
      type = "binary",
      notes = "Screened (Methods, Covariate screening) but not retained in the final model (Eqs. 1-2)."
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened (Table 1: median 19.7, range 9.5-88 kg) but not retained in the final model."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Screened (Table 1: median 18.3, range 7.6-46.8 g/L) but not retained."
    ),
    HCT = list(
      description = "Haematocrit",
      units = "%",
      type = "continuous",
      notes = "Screened (Table 1: median 41.2 %; the printed upper range 447.6 is a misprint) but not retained."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 139L,
    n_studies = 1L,
    n_observations = 432L,
    age_range = "1.1-15.6 years",
    age_median = "5.3 years",
    weight_range = "9.5-88 kg",
    weight_median = "19.7 kg",
    sex_female_pct = 73.4,
    race_ethnicity = "Chinese (single centre, Guangzhou Women and Children's Medical Center).",
    disease_state = paste0(
      "Children with steroid-resistant or steroid-dependent (refractory) ",
      "nephrotic syndrome diagnosed per the 2012 KDIGO guideline, onset ",
      "before 18 years, on oral tacrolimus (Prograf) plus prednisone or ",
      "prednisolone, August 2013 to December 2018."
    ),
    dose_range = paste0(
      "Starting dose 500-3000 ug orally every 12 h, titrated to a whole-blood ",
      "trough of 5-10 ug/L."
    ),
    regions = "China (Guangzhou).",
    co_medication = "Steroids in all patients; Wuzhi capsule 6/139 (Table 1).",
    genotypes = paste0(
      "CYP3A5 rs776746 *1/*1 18, *1/*3 58, *3/*3 63; CTLA4 rs4553808 GG 2, ",
      "GA 29, AA 108 (Table 2). Ten further podocyte / immune SNPs were ",
      "screened and not retained."
    ),
    assay = "Whole-blood tacrolimus by EMIT (Viva-E, Siemens).",
    notes = paste0(
      "432 concentrations, of which 35 are peak samples (0.5-3 h post-dose) ",
      "and the rest pre-dose troughs. Table 1 prints gender as 37/102 ",
      "(male/female). Phoenix NLME 8.1, FOCE-ELS; 1000-replicate bootstrap."
    )
  )

  ini({
    # Huang 2022 Table 3 'Parameter estimation results and bootstrap results
    # of the final population pharmacokinetic model'; final-model equations
    # printed in Results as Eqs. 1-2:
    #   V (L) = 192.03 x exp(etaV)
    #   CL (L/h) = 10.54 x (AGE/5.3)^0.31 x exp(-0.34 x CTLA4 GA)
    #     x exp(-0.15 x CTLA4 AA) x exp(-0.25 x CYP3A5*3/*3)
    #     x exp(-0.34 x combWZ) x exp(etaCL)
    lka <- fixed(log(4.48)); label("Absorption rate constant ka (1/h)") # Methods 'Base model': Ka fixed at 4.48 1/h from the literature
    lcl <- log(10.54); label("Apparent oral clearance CL/F at 5.3 years, CTLA4 GG, CYP3A5 expresser, no Wuzhi capsule (L/h)") # Table 3 tvCL = 10.54 L/h; Eq. 2
    lvc <- log(192.03); label("Apparent volume of distribution V/F (L)") # Table 3 tvV = 192.03 L; Eq. 1

    e_age_cl <- 0.31; label("Power exponent of age on CL/F, reference 5.3 years (unitless)") # Table 3 dCLdAGE = 0.31; Eq. 2
    e_ctla4_ga_cl <- -0.34; label("Exponential effect of CTLA4 rs4553808 GA (vs GG) on CL/F (unitless)") # Table 3 dCLdCTLA4 GA = -0.34; Eq. 2
    e_ctla4_aa_cl <- -0.15; label("Exponential effect of CTLA4 rs4553808 AA (vs GG) on CL/F (unitless)") # Table 3 dCldCTLA4 AA = -0.15; Eq. 2
    e_cyp3a5_expr_cl <- -0.25; label("Exponential effect of CYP3A5*3/*3 non-expresser status on CL/F, applied to (1 - CYP3A5_EXPR) (unitless)") # Table 3 dCLdCYP3A5 *3/*3 = -0.25; Eq. 2
    e_conmed_wuzhi_cl <- -0.34; label("Exponential effect of concomitant Wuzhi capsule on CL/F (unitless)") # Table 3 dCLdCombWZ = -0.34; Eq. 2

    # Table 3 rows 'omega2 V' = 1.13 (shrinkage 38.09%) and 'omega2 CL' = 0.13
    # (shrinkage 30.02%), labelled in the footnote as variances of the
    # interindividual variability; no covariance is reported.
    etalvc ~ 1.13 # Table 3 omega2 V = 1.13 (variance)
    etalcl ~ 0.13 # Table 3 omega2 CL = 0.13 (variance)

    # Results: 'proportional models were chosen to describe the ... residual
    # variability'; Methods Yobs = Ypred x (1 + eps1).
    propSd <- 0.2995; label("Proportional residual error (fraction)") # Table 3 sigma = 29.95% (bootstrap median 0.30)
  })

  model({
    ka <- exp(lka)
    # Eq. 2. The CYP3A5 coefficient multiplies the source *3/*3 indicator,
    # which is 1 - CYP3A5_EXPR.
    cl <- exp(lcl + etalcl) * (AGE / 5.3)^e_age_cl *
      exp(e_ctla4_ga_cl * SNP_CTLA4_RS4553808_GA) *
      exp(e_ctla4_aa_cl * SNP_CTLA4_RS4553808_AA) *
      exp(e_cyp3a5_expr_cl * (1 - CYP3A5_EXPR)) *
      exp(e_conmed_wuzhi_cl * CONMED_WUZHI)
    # Eq. 1
    vc <- exp(lvc + etalvc)
    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Dose in mg, vc in L -> mg/L; x 1000 gives ng/mL (= ug/L)
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
