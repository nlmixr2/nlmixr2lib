Li_2017_CC292 <- function() {
  description <- paste(
    "Two-compartment population PK model for oral CC-292 (spebrutinib, a",
    "potent Bruton tyrosine kinase inhibitor) in 145 pooled subjects: 32",
    "healthy adults (AVL-292-004) and 113 patients with relapsed and/or",
    "refractory B-cell malignancies including chronic lymphocytic leukemia",
    "(AVL-292-003). First-order absorption with a single absorption lag,",
    "linear elimination from the central compartment, with linear-deviation",
    "female-sex effect on apparent clearance (females have 26% lower CL/F)",
    "and a power age effect on apparent central volume (reference age 62",
    "years). Residual variability is split into healthy-volunteer and",
    "patient strata."
  )
  reference <- paste(
    "Li Y, Ramirez-Valle F, Xue Y, Ventura JI, Gouedard O, Mei J,",
    "Takeshita K, Palmisano M, Zhou S. Population Pharmacokinetics and",
    "Exposure Response Assessment of CC-292, a Potent BTK Inhibitor, in",
    "Patients With Chronic Lymphocytic Leukemia. J Clin Pharmacol.",
    "2017;57(10):1279-1289. doi:10.1002/jcph.923",
    sep = " "
  )
  vignette <- "Li_2017_CC292"
  units    <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "Linear-deviation categorical effect on CL/F: multiplier",
        "(1 - 0.26 * SEXF) per Li 2017 final covariate model (page 1283,",
        "covariate equations following Table 2). Females have ~26% lower",
        "apparent clearance than males; the paper concluded this effect is",
        "not clinically relevant. Reference cohort sex split: 60% male,",
        "40% female (Li 2017 Table 1)."
      ),
      source_name        = "SEX"
    ),
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power covariate on apparent central volume V2/F:",
        "V2/F = 158 * (AGE / 62)^0.946 per Li 2017 final covariate",
        "equations (page 1283, following Table 2). Reference 62 years is",
        "the cohort median age (Li 2017 Table 1). Time-fixed at baseline.",
        "The paper concluded the age effect on V2/F is not clinically",
        "relevant. Cohort age range 20-89 years."
      ),
      source_name        = "AGE"
    ),
    DIS_HEALTHY = list(
      description        = "Healthy-participant cohort indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (patient with B-cell malignancy)",
      notes              = paste(
        "Stratifies the residual error: HNP (healthy normal participants",
        "from AVL-292-004) have a lower residual variance (sigma^2 = 0.234)",
        "than patients with relapsed/refractory B-cell malignancies from",
        "AVL-292-003 (sigma^2 = 0.659), reflecting the better data quality",
        "of the well-controlled healthy-adult study. Cohort split: 32",
        "healthy participants (22.1%) and 113 patients (77.9%) of the 145",
        "total subjects (Li 2017 Methods 'Patients and Study Design'",
        "section, study counts; Table 1 demographics)."
      ),
      source_name        = "HNP"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 145L,
    n_studies      = 2L,
    n_observations = 3156L,
    age_range      = "20-89 years (median 62)",
    age_median     = "62 years",
    weight_range   = "49.9-128.4 kg (median 79.5)",
    weight_median  = "79.5 kg",
    height_range   = "149-200 cm (median 170)",
    bmi_range      = "19-41 kg/m^2 (median 27)",
    sex_female_pct = 40,
    race_ethnicity = c(White = 83.4, Black = 13.8, Asian = 1.4, Other = 1.4),
    ethnicity      = c("Hispanic or Latino" = 10.3, "Not Hispanic or Latino" = 89.7),
    disease_state  = paste(
      "Pooled cohort: 32 healthy adults from AVL-292-004 (Part 1, 7-day",
      "QD multiple-dose) and 113 patients from AVL-292-003 with relapsed",
      "and/or refractory B-cell malignancies (B-NHL, CLL/SLL,",
      "Waldenstrom's macroglobulinemia). 25/145 subjects (17.2%) had",
      "moderate renal impairment (CLcr 30-60 mL/min); none had severe",
      "renal impairment."
    ),
    dose_range     = paste(
      "Oral solid dosage. Healthy adults: 50, 100, or 200 mg QD x 7 days",
      "(AVL-292-004 Part 1). Patients: QD doses of 125, 250, 400, 625,",
      "750, 1000 mg or BID doses of 375 or 500 mg administered in 28-day",
      "cycles (AVL-292-003)."
    ),
    regions        = "Not stated; both studies sponsored by Avila Therapeutics / Celgene Corporation",
    cohort_split   = "32 healthy (22.1%) + 113 patients (77.9%) = 145 total",
    renal_impair_pct = 17.2,
    notes          = paste(
      "Baseline demographics from Li 2017 Table 1 (pooled column). The",
      "non-white subjects were grouped together for the race covariate",
      "analysis (race was not retained in the final model).",
      "The bioanalytical LLOQ was 0.5 ng/mL; the model characterised",
      "concentrations in the log range -0.67 to 9.152, i.e. 0.51 ng/mL",
      "to 9440 ng/mL (Li 2017 Results, Model Evaluation)."
    )
  )

  ini({
    # Structural PK parameters from Li 2017 Table 2 (final model column).
    # Typical values correspond to the reference subject: male (SEXF = 0)
    # with age = 62 years (the cohort median age).
    lka     <- log(0.974); label("First-order absorption rate constant (ka, 1/h)")                     # Li 2017 Table 2 (ka = 0.974 1/h)
    lcl     <- log(134);   label("Apparent clearance for reference male subject (CL/F, L/h)")          # Li 2017 Table 2 (CL/F = 134 L/h)
    lvc     <- log(158);   label("Apparent central volume of distribution at AGE = 62 years (V2/F, L)") # Li 2017 Table 2 (V2/F = 158 L)
    lq      <- log(18.7);  label("Apparent inter-compartmental clearance (Q/F, L/h)")                  # Li 2017 Table 2 (Q/F = 18.7 L/h)
    lvp     <- log(72);    label("Apparent peripheral volume of distribution (V3/F, L)")               # Li 2017 Table 2 (V3/F = 72 L)
    lalag   <- log(0.427); label("Absorption lag time (Alag1, h)")                                     # Li 2017 Table 2 (Alag1 = 0.427 h)

    # Covariate effects (Li 2017 final covariate model equations, page 1283).
    # Sex effect on CL/F: categorical linear-deviation form per Methods
    # ("Categorical covariates" formula). Females (SEXF = 1) have 26% lower
    # CL/F: CL/F = 134 * (1 - 0.26 * SEXF). Coefficient sign is negative to
    # encode the reduction; the paper's Table 2 entry "Effect of sex on
    # CL/F = 0.26" reports the magnitude.
    e_sexf_cl <- -0.26;  label("Linear-deviation effect of female sex on CL/F (fractional)")           # Li 2017 Table 2 (Effect of sex on CL/F = 0.26 magnitude; females -26%)

    # Age effect on V2/F: power form (preferred over linear by OFV decrease
    # -2758 -> -2722 in stepwise covariate selection). Reference age 62
    # years (cohort median). V2/F = 158 * (AGE / 62)^0.946.
    e_age_vc  <-  0.946; label("Power exponent of (AGE / 62) on V2/F (unitless)")                      # Li 2017 Table 2 (Effect of age on V2/F = 0.946)

    # IIV (Li 2017 Table 2: variance/covariance reported on log-eta scale).
    #   omega^2 V2/F           = 0.755                    -> CV = sqrt(exp(0.755) - 1) = 87%
    #   omega^2 CL/F           = 0.317                    -> CV = sqrt(exp(0.317) - 1) = 60%
    #   omega(V2/F):omega(CL/F) = 0.328 (covariance)
    #   correlation rho        = 0.328 / sqrt(0.317 * 0.755) = 0.670
    # The paper reports IIV magnitudes in the Discussion as ~56% on CL/F and
    # ~86% on V2/F (final model); the variance-to-CV check is consistent.
    # nlmixr2 block order: etalcl + etalvc -> c(var_cl, cov, var_vc).
    etalcl + etalvc ~ c(0.317, 0.328, 0.755)                                                          # Li 2017 Table 2 (IIV block: omega^2 CL/F, covariance, omega^2 V2/F)

    # Residual error. The source data were ln-transformed and the residual
    # was modelled as additive on the log scale with variance sigma^2; this
    # is equivalent to a proportional residual model on the linear scale
    # where propSd = sqrt(sigma^2). Variance is stratified by study
    # population (Li 2017 Methods, Population Pharmacokinetic Model Building).
    propSd_hv <- sqrt(0.234); label("Proportional residual SD for healthy participants (fraction)")    # Li 2017 Table 2 (HNP sigma^2 = 0.234)
    propSd_pt <- sqrt(0.659); label("Proportional residual SD for patients (fraction)")                # Li 2017 Table 2 (Patients sigma^2 = 0.659)
  })

  model({
    # Individual PK parameters. CL/F is shifted by the SEXF linear-deviation
    # effect; V2/F is scaled by the power AGE/62 factor. Both reference
    # categories (male sex, age 62 years) collapse to the published typical
    # values 134 L/h and 158 L (Li 2017 Discussion: "Typical values of
    # CC-292 CL/F and V2/F for male subjects with a median age of 62 years
    # were 134 L/h and 158 L").
    ka     <- exp(lka)
    cl     <- exp(lcl + etalcl) * (1 + e_sexf_cl * SEXF)
    vc     <- exp(lvc + etalvc) * (AGE / 62)^e_age_vc
    q      <- exp(lq)
    vp     <- exp(lvp)
    alag_d <- exp(lalag)

    # Micro-constants for the 2-compartment central-disposition model.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ODE system.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Absorption lag.
    alag(depot) <- alag_d

    # Plasma CC-292 concentration. Dose in mg, vc in L -> central / vc in
    # mg/L = ug/mL; rescaled by 1000 to ng/mL to match Li 2017 reporting
    # (300 ng/mL threshold in the exposure-response analysis; 0.5 ng/mL
    # bioanalytical LLOQ).
    Cc <- 1000 * central / vc

    # Population-stratified proportional residual error. HNP cohort uses
    # propSd_hv (sigma^2 = 0.234) and patient cohort uses propSd_pt
    # (sigma^2 = 0.659) per Li 2017 Table 2 (Residual variability rows).
    propSdEff <- propSd_hv * DIS_HEALTHY + propSd_pt * (1 - DIS_HEALTHY)
    Cc ~ prop(propSdEff)
  })
}
