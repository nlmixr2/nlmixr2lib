Charles_2007_tafenoquine <- function() {
  description <- "One-compartment first-order-absorption population PK model for oral tafenoquine in adult Australian soldiers on weekly malaria prophylaxis (Charles 2007)"
  reference <- "Charles BG, Miller AK, Nasveld PE, Reid MG, Harris IE, Edstein MD. Population pharmacokinetics of tafenoquine during malaria prophylaxis in healthy subjects. Antimicrob Agents Chemother. 2007;51(8):2709-2715. doi:10.1128/AAC.01183-06"
  vignette <- "Charles_2007_tafenoquine"
  units <- list(time = "hour", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed (study baseline). Enters CL/F and V/F as a centered linear effect, parameterized as (1 + theta * WT/80.9) with the cohort-mean weight of 80.9 kg as the centering constant (Charles 2007 Table 1 footnote d, Table 2 footnote a).",
      source_name        = "WT"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age in years",
      units       = "years",
      type        = "continuous",
      notes       = "Screened as a centered linear effect (AGE/25.4) on CL/F and V/F. Age on V/F was significant (delta-OFV = -9) but positively correlated with weight and therefore not retained in the final model (Charles 2007 Table 1 models 2-3; Results)."
    ),
    CRCL = list(
      description = "Estimated creatinine clearance (Cockcroft-Gault)",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Screened as a centered linear effect (CLCR/121) on CL/F; not significant (delta-OFV = -4) and not retained (Charles 2007 Table 1 model 4)."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened on CL/F (delta-OFV = -3) and V/F (delta-OFV = -12) but not retained in the final model (Charles 2007 Table 1 models 7-8; Results). Cohort was 476 male / 14 female."
    ),
    PHOS = list(
      description = "Phospholipidosis-present indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened in the 77-subject phospholipidosis substudy; not significant on CL/F or V/F (Charles 2007 Table 1 models 5-6). Not a registered canonical covariate; recorded here for source-trace completeness only."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 490,
    n_studies      = 1,
    age_range      = "18-47 years",
    age_median     = "25.4 years (mean +/- 5.3 SD)",
    weight_range   = "50-135 kg",
    weight_median  = "80.9 kg (mean +/- 11.9 SD)",
    sex_female_pct = 14 / 490 * 100,
    race_ethnicity = c(White = (490 - 8) / 490 * 100, Other = 8 / 490 * 100),
    disease_state  = "Healthy adult Australian soldiers on malaria prophylaxis during 6-month deployment to East Timor",
    dose_range     = "200 mg PO once daily for 3 days (loading) then 200 mg PO once weekly for ~6 months",
    regions        = "Australia (deployed to East Timor)",
    notes          = "Phase III prospective randomised double-blind trial of weekly tafenoquine for malaria prophylaxis. Baseline demographics in Charles 2007 Results paragraph 1: 476 males and 14 females; 482 / 490 subjects of Caucasian background. 1,925 plasma concentration-time points contributed to the popPK fit. Glucose-6-phosphate dehydrogenase-normal subjects only."
  )

  ini({
    # Structural parameters. Final-model partials (theta_1..theta_5)
    # are reported in Charles 2007 Table 2. The typical CL/F and V/F at the
    # cohort-mean weight of 80.9 kg are theta_1*(1+theta_4) = 4.37 L/h and
    # theta_2*(1+theta_5) = 1901 L respectively, matching the Results text.
    lka     <- log(0.243);  label("Absorption rate constant (Ka, 1/h)")            # Charles 2007 Table 2 (theta_3 = 0.243 1/h)
    lcl     <- log(3.02);   label("Clearance partial theta_1 (CL/F, L/h)")          # Charles 2007 Table 2 (theta_1 = 3.02 L/h; typical CL/F at 80.9 kg = 4.37 L/h)
    lvc     <- log(1110);   label("Central volume partial theta_2 (V/F, L)")        # Charles 2007 Table 2 (theta_2 = 1110 L; typical V/F at 80.9 kg = 1901 L)

    # Centered linear weight effects, parameterized as (1 + theta * WT/80.9).
    # See Charles 2007 Table 2 footnote a and Results paragraph 3.
    e_wt_cl <- 0.448;       label("Linear coefficient of WT/80.9 on CL/F (unitless)")  # Charles 2007 Table 2 (theta_4 = 0.448)
    e_wt_vc <- 0.713;       label("Linear coefficient of WT/80.9 on V/F (unitless)")   # Charles 2007 Table 2 (theta_5 = 0.713)

    # IIV. omega^2 = log(CV^2 + 1):
    #   18% CV -> log(1 + 0.18^2) = 0.03189
    #   22% CV -> log(1 + 0.22^2) = 0.04727
    #   76% CV -> log(1 + 0.76^2) = 0.45603
    # The paper estimated a covariance between eta_CL/F and eta_V/F (Results
    # paragraph 4: OFV dropped 22,265 -> 22,248 vs the diagonal-only model)
    # but did not report the numeric value, so independent diagonals are used
    # here; see vignette Assumptions and deviations.
    etalcl ~ 0.03189   # Charles 2007 Table 2 (IIV CL/F = 18% CV)
    etalvc ~ 0.04727   # Charles 2007 Table 2 (IIV V/F = 22% CV)
    etalka ~ 0.45603   # Charles 2007 Table 2 (IIV Ka = 76% CV)

    # Residual error: combined proportional + additive (Charles 2007 Methods
    # "Population pharmacokinetic modeling" paragraph 4; values in Table 2).
    # Additive 22.9 ng/mL = 0.0229 ug/mL = 0.0229 mg/L.
    propSd <- 0.059;        label("Proportional residual error (fraction)")            # Charles 2007 Table 2 (RUV 5.9% CV)
    addSd  <- 0.0229;       label("Additive residual error (ug/mL)")                   # Charles 2007 Table 2 (RUV 22.9 ng/mL)
  })

  model({
    # Centered linear weight effects (cohort mean = 80.9 kg, Table 2 footnote a)
    cl <- exp(lcl + etalcl) * (1 + e_wt_cl * WT / 80.9)
    vc <- exp(lvc + etalvc) * (1 + e_wt_vc * WT / 80.9)
    ka <- exp(lka + etalka)

    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Plasma concentration: dose in mg, volume in L -> mg/L = ug/mL.
    # Paper reports concentrations in ng/mL (1 ug/mL = 1000 ng/mL).
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
