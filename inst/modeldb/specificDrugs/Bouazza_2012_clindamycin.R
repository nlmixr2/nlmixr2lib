Bouazza_2012_clindamycin <- function() {
  description <- paste(
    "One-compartment population PK model for clindamycin administered",
    "orally (immediate-release tablet) or intravenously (20-min infusion)",
    "in 50 adult patients (ages 18-93 y, body weight 23-133 kg) treated",
    "for bone and joint infections (Bouazza 2012). First-order absorption",
    "for oral dosing with estimated absolute bioavailability F = 0.876;",
    "apparent clearance CL/F = 15.2 L/h at 70 kg, with an estimated",
    "(non-allometric) body-weight exponent of 0.497 on CL. Apparent",
    "volume V/F = 66.2 L and absorption rate Ka = 0.967 1/h, neither",
    "carrying retained interindividual variability. IIV is retained only",
    "on CL/F (omega = 0.39). Residual variability is proportional (sigma",
    "= 0.38). Rifampicin co-administration was screened but not retained",
    "in the final model (see covariatesDataExcluded)."
  )
  reference <- paste(
    "Bouazza N, Pestre V, Jullien V, Curis E, Urien S, Salmon D,",
    "Treluyer JM.",
    "Population pharmacokinetics of clindamycin orally and intravenously",
    "administered in patients with osteomyelitis.",
    "British Journal of Clinical Pharmacology. 2012;74(6):971-977.",
    "doi:10.1111/j.1365-2125.2012.04292.x.",
    sep = " "
  )
  vignette <- "Bouazza_2012_clindamycin"
  units    <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Used as a continuous covariate on apparent clearance via the",
        "estimated power model CL/F = theta_CL * (BW/70)^cov_WT (Bouazza",
        "2012 Results, 'Population pharmacokinetics'). Cohort range 23",
        "to 133 kg (mean 70 kg, 95% CI 45-113 kg). The reference",
        "weight is 70 kg (the typical adult)."
      ),
      source_name        = "BW"
    )
  )

  covariatesDataExcluded <- list(
    RIFAMPICIN = list(
      description        = "Rifampicin co-administration (binary)",
      units              = "(binary)",
      type               = "binary",
      reference_category = NULL,
      notes              = paste(
        "Co-treatment with rifampicin was tested as a categorical effect",
        "on clindamycin clearance (4 of 50 patients co-treated, 8%).",
        "The effect (43% increase in CL/F, SE 0.17) was significant by",
        "the Wald test (P = 0.04) but not by the more conservative",
        "likelihood-ratio test, and the small co-treated cohort (n = 4)",
        "did not provide enough power to retain it. The paper therefore",
        "removed rifampicin from the final model (Bouazza 2012 Results,",
        "'Population pharmacokinetics'; Discussion notes the mechanism",
        "as CYP3A4 induction by rifampicin)."
      )
    ),
    AGE = list(
      description        = "Age",
      units              = "year",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Tabulated in the cohort summary (Table 1: 56.7 +/- 3.0 years,",
        "range 18-93). Bouazza 2012 Methods 'Modelling strategy and",
        "population pharmacokinetic model' lists age among the screened",
        "covariates; only body weight was retained in the final model."
      )
    ),
    RENALFAIL = list(
      description        = "Renal failure (binary; Cockcroft-Gault estimated CrCL)",
      units              = "(binary)",
      type               = "binary",
      reference_category = NULL,
      notes              = paste(
        "Tabulated in the cohort summary (Table 1: 5 of 50 patients,",
        "10%). Renal failure was assessed because clindamycin dosage is",
        "not adjusted in renal insufficiency per the label; no",
        "significant effect on clindamycin clearance was retained in",
        "the final model (Bouazza 2012 Discussion)."
      )
    ),
    HEPATICFAIL = list(
      description        = "Hepatic failure (binary)",
      units              = "(binary)",
      type               = "binary",
      reference_category = NULL,
      notes              = paste(
        "Tabulated in the cohort summary (Table 1: 2 of 50 patients,",
        "4.2%). Hepatic function was screened because clindamycin is",
        "primarily hepatically metabolised; no significant effect was",
        "retained in the final model."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 50L,
    n_studies      = 1L,
    n_observations = 122L,
    age_range      = "18-93 years (Table 1 mean 56.7)",
    weight_range   = "23-133 kg (Table 1 mean 69.9; Methods reports a 95% reference range of 45-113 kg)",
    sex_female_pct = 40,
    race_ethnicity = NA,
    disease_state  = "Adult patients with bone and joint infections (osteomyelitis)",
    dose_range     = paste(
      "600 mg three times daily (q8h) by oral tablet or 20-min",
      "intravenous infusion; one patient received 600 mg four times",
      "daily and two patients received 600 mg once daily. Median",
      "weight-normalised daily dose 26.0 mg/kg/day (range 13.5-41.8).",
      "Of 50 patients, 22 received only oral, 22 only i.v., and 6",
      "received both routes (58 PO plus 64 i.v. plasma concentrations).",
      "All patients sampled at steady state."
    ),
    regions        = "France (Cochin Hospital, Paris; retrospective therapeutic drug monitoring, 2008-2010)",
    notes          = paste(
      "Two of 122 plasma concentrations were below the 0.1 mg/L limit",
      "of quantification and were handled as left-censored data by the",
      "Monolix M3 method. BMI distribution: 24% normal weight (< 25",
      "kg/m^2), 68% overweight (25-30 kg/m^2), and 8% obese (> 30",
      "kg/m^2). Co-medications include rifampicin in 4 patients (8%);",
      "see covariatesDataExcluded. Race, ethnicity, and other demographic",
      "details are not reported in the source. Model fit with Monolix",
      "version 4 (SAEM algorithm + MCMC; 10 MCMC chains). Validation",
      "via prediction-corrected VPC and NPDE."
    )
  )

  ini({
    # Structural parameters - typical values at the reference body weight of 70 kg.
    # Bouazza 2012 Table 2 ('Structural model' block); apparent values are
    # CL/F and V/F derived from i.v. infusion data so that F = 1 on the i.v.
    # side and F = 0.876 on the oral side.
    lka <- log(0.967); label("Absorption rate constant (Ka, 1/h)")                      # Table 2: Ka = 0.967 1/h (RSE 26%)
    lcl <- log(15.2);  label("Clearance at WT = 70 kg (CL, L/h)")                       # Table 2: CL = 15.2 L/h (RSE 8%)
    lvc <- log(66.2);  label("Volume of distribution (V, L)")                           # Table 2: V = 66.2 L (RSE 9%)

    # Estimated body-weight exponent on CL (NOT canonical 0.75 allometric).
    # Bouazza 2012 Results 'Population pharmacokinetics' final covariate model:
    #   CL = theta_CL * (BW / 70)^cov_WT
    # Table 2 reports cov_WT = 0.497 (RSE 36%); the paper estimated this
    # exponent (it is not a fixed allometric power) and notes the addition of
    # the BW covariate "decreased the AIC/BIC criteria, resulted in a 5.86
    # units decrease in the objective function value and improved the
    # goodness of fit".
    e_wt_cl <- 0.497; label("Estimated body-weight exponent on CL (unitless)")          # Table 2: cov_WT = 0.497 (RSE 36%)

    # Bioavailability of the oral form (estimated using a logit-normal
    # transformation for boundedness; reported as F = 0.876 with RSE 11%).
    # Logit-normal here refers to the parameterisation used during
    # estimation, not to IIV: no IIV on F was retained (Bouazza 2012
    # Results: "BSVs were described by an exponential error model and
    # retained only for apparent clearance"). Encoded with the canonical
    # log-transformed depot-fraction parameter so f(depot) = exp(lfdepot)
    # = 0.876 in the model() block.
    lfdepot <- log(0.876); label("Oral bioavailability F (unitless; applied to depot only)") # Table 2: F = 0.876 (RSE 11%)

    # IIV - exponential model q_i = q_pop * exp(eta_i), eta_i ~ N(0, omega^2).
    # Only IIV on apparent clearance was retained (Bouazza 2012 Results).
    # Monolix reports omega as the SD of eta on the log scale; nlmixr2's ~
    # syntax expects the variance, so we encode 0.39^2 = 0.1521.
    etalcl ~ 0.1521  # Table 2: omega_CL/F = 0.39 (RSE 12%); variance = 0.39^2 = 0.1521

    # Residual error - proportional model selected over proportional + additive
    # by likelihood-ratio testing (Bouazza 2012 Methods 'Modelling strategy').
    # Monolix reports sigma as the proportional SD on the linear scale, which
    # maps directly onto nlmixr2's propSd.
    propSd <- 0.38; label("Proportional residual error (fraction)")                    # Table 2: sigma = 0.38 (RSE 8%)
  })

  model({
    # 1. Derived covariate terms
    # Apparent clearance scales with body weight via the estimated power
    # exponent e_wt_cl around the reference body weight of 70 kg. The
    # exponent is 0.497 (Table 2), which is closer to the 0.5 reported
    # for some adult cohorts than to the canonical theoretical 0.75.

    # 2. Individual PK parameters
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc)

    # 3. Micro-constants
    kel <- cl / vc

    # 4. ODE system: depot -> central, first-order absorption, linear elimination
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # 5. Bioavailability - F applies to the oral (depot) input only.
    # Intravenous doses bypass the depot compartment by going directly
    # to central, so the F = 0.876 attenuation is only experienced by
    # oral doses.
    f(depot) <- exp(lfdepot)

    # 6. Observation: plasma clindamycin concentration (dose in mg, V in L
    # -> mg/L, matching the source's reporting units).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
