Hu_2017_peginterferon_beta_1a <- function() {
  description <- "One-compartment population PK model for peginterferon beta-1a in adults with relapsing multiple sclerosis (Hu 2017). First-order SC absorption with the absorption rate constrained above the elimination rate to avoid flip-flop kinetics. BMI is a covariate on both clearance and volume of distribution."
  reference <- "Hu X, Hang Y, Cui Y, Zhang J, Liu S, Seddighzadeh A, Deykin A, Nestorov I. Population-Based Pharmacokinetic and Exposure-Efficacy Analyses of Peginterferon Beta-1a in Patients With Relapsing Multiple Sclerosis. J Clin Pharmacol. 2017;57(8):1005-1016. doi:10.1002/jcph.883"
  vignette <- "Hu_2017_peginterferon_beta_1a"
  units <- list(time = "hour", dosing = "ug", concentration = "ng/mL")

  covariateData <- list(
    BMI = list(
      description        = "Body mass index at baseline",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at baseline. Used with a power form on clearance ((BMI / 23.71)^e_bmi_cl) and an exponential-deviation form on central volume (exp(e_bmi_vc * (BMI - 23.71))). Reference value 23.71 kg/m^2 is the typical BMI used in Hu 2017 equations 12-13.",
      source_name        = "BMI"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 809,
    n_studies      = 1,
    age_range      = "20.5-54.7 years (2.5th-97.5th percentile; Table 1)",
    age_median     = "36.6 years",
    weight_range   = "46.0-103 kg (2.5th-97.5th percentile; Table 1)",
    weight_median  = "65.0 kg",
    bmi_range      = "17.4-35.5 kg/m^2 (2.5th-97.5th percentile; Table 1)",
    bmi_median     = "23.3 kg/m^2",
    sex_female_pct = 70.5,
    race_ethnicity = c(White = 82.6, Asian = 11.9, Other = 5.6),
    disease_state  = "Relapsing multiple sclerosis (ADVANCE phase 3 trial).",
    dose_range     = "125 ug SC every 2 weeks or every 4 weeks",
    regions        = "26 countries, 183 sites (multicenter)",
    trial          = "NCT00906399; ATTAIN extension NCT00910689 referenced in Hu 2017",
    notes          = "Hu 2017 Table 1 final PK population: 239 males and 570 females, predominantly White (n=668) and Asian (n=96). 1512 patients were randomized 1:1:1 to placebo, peginterferon beta-1a 125 ug SC Q2W, or peginterferon beta-1a 125 ug SC Q4W in ADVANCE; the final PK analysis excluded BLQ data (62%), concentrations beyond 10 days postdose due to missing dose information (4%), three sparse PK subjects with positive baseline measurement, one outlier concentration (14700 pg/mL), and concentrations with positive anti-IFN antibodies (0.8%). Intensive PK sampling was performed in 25 subjects (12 on Q2W, 13 on Q4W); the remainder had sparse sampling."
  )

  ini({
    # Structural parameters at typical BMI = 23.71 kg/m^2 (Hu 2017 Table 3 final model)
    lcl         <- log(3.28);  label("Typical clearance at BMI 23.71 kg/m^2 (L/h)")          # Hu 2017 Table 3 theta_1
    lvc         <- log(435);   label("Typical central volume at BMI 23.71 kg/m^2 (L)")       # Hu 2017 Table 3 theta_2
    ltheta_diff <- log(0.207); label("Log of (ka - kel) used to constrain absorption rate above elimination rate (1/h)") # Hu 2017 Table 3 theta_5 (eq 1)

    # BMI covariate effects (final model only retains BMI; Hu 2017 eqs 12-13)
    e_bmi_cl <- 0.779;  label("Exponent of BMI/23.71 on clearance (unitless)")                # Hu 2017 Table 3 theta_3
    e_bmi_vc <- 0.0353; label("Coefficient of (BMI - 23.71) inside exp() on central volume (per kg/m^2)") # Hu 2017 Table 3 theta_4

    # Inter-individual variability (omega^2, log-normal IIV on CL and Vc)
    etalcl ~ 0.145  # Hu 2017 Table 3 omega^2_CL (~40 percent CV, 63 percent shrinkage)
    etalvc ~ 0.352  # Hu 2017 Table 3 omega^2_V  (~65 percent CV, 57 percent shrinkage)

    # Residual error: NONMEM log-additive form Y = LOG(F) + SD1*EPS with $SIGMA 1 FIXED,
    # i.e. log-normal residual with SD on the log-scale = SD1.
    expSd <- 0.566; label("Log-scale residual SD (unitless)")                                  # Hu 2017 Table 3 SD1
  })

  model({
    # Individual PK parameters with BMI covariate effects (Hu 2017 eqs 12-13)
    cl <- exp(lcl + etalcl) * (BMI / 23.71)^e_bmi_cl
    vc <- exp(lvc + etalvc) * exp(e_bmi_vc * (BMI - 23.71))

    # Absorption rate constrained above elimination rate (Hu 2017 eq 1):
    #   ka_i = theta_diff + cl_i / vc_i, with theta_diff > 0 to prevent flip-flop.
    kel <- cl / vc
    ka  <- exp(ltheta_diff) + kel

    # Bioavailability is fixed at F = 1 because no IV data were available
    # (Hu 2017 Methods, Population PK Model section). No explicit lfdepot parameter
    # is needed since the rxode2 default for depot is F = 1.

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Concentration: dose in ug, vc in L -> ug/L = ng/mL
    Cc <- central / vc
    Cc ~ lnorm(expSd)
  })
}
